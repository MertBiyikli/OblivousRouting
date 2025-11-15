//
// Created by Mert Biyikli on 19.09.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_FRAMEWORK_H
#define OBLIVIOUSROUTING_RAECKE_FRAMEWORK_H


#include <memory>
#include <vector>
#include <unordered_map>
#include "../solver/solver.h"



/**
 * This class implements Räcke’s oblivious‐routing solver
 * as in his STOC 2008 paper using the FRT algorithm.
 * The templates are:
 * - Algorithm: the RaeckeBase algorithm (e.g., RaeckeFRT or RaeckeMST)
 * - Transform: the RaeckeTransformBase (e.g., RaeckeFRTTransform or RaeckeMSTTransform)
 */

template<typename Algorithm, typename Transform>
class RaeckeFramework : public ObliviousRoutingSolver{
    Algorithm m_algorithm;
    Transform m_transform;

public:
    void solve(const Graph & g) override{
        m_algorithm.setGraph(g);
        m_algorithm.run();
        // scale the flow to meet unit flow

        iteration_count = m_algorithm.getTrees().size();
        this->oracle_running_times = m_algorithm.oracle_running_times;
        this->pure_oracle_running_times = m_algorithm.pure_oracle_running_times;
    }

    void storeFlow() override {

        for (size_t i = 0; i < m_algorithm.getTrees().size(); ++i) {
            m_transform.addTree(m_algorithm.getTrees()[i],  m_algorithm.getLambdas()[i], m_algorithm.getGraphs()[i]);
        }



        // store the flow
        // given the demand map
        auto const& routingRaecke = m_transform.getRoutingTable();
        for (const auto& [edge, demandMap] : routingRaecke) {
            for (const auto& [d, fraction] : demandMap) {
                if (d.first > d.second) continue; // skip trivial cases
                f_e_st[edge][d]=fraction;
            }
        }

        scaleDownFlow();
    }

    double getCongestion() {
        double max_congestion = 0.0;

        auto const& routingRaecke = m_transform.getRoutingTable();

        for (const auto& [edge, demandMap] : routingRaecke) {
            double total_flow = 0.0;
            for (const auto& [d, fraction] : demandMap) {
                total_flow += fraction;
            }

            double capacity = m_algorithm.getGraph().getEdgeCapacity(edge.first, edge.second);
            if (capacity > 0) {
                double congestion = total_flow / capacity;
                if (congestion > max_congestion) {
                    max_congestion = congestion;
                }
            }
        }

        return max_congestion;
    }

    void scaleDownFlow() {
        // scale the flow to meet unit flow
        std::unordered_map<std::pair<int, int>, double > outgoingflow_per_commodity;

        for ( const auto& [edge, flowMap]:f_e_st) {
            for (const auto& [com, flow_value]  : flowMap) {
                if ( flow_value < 1e-15 ) continue; // ignore zero flows
                if (!outgoingflow_per_commodity.contains(com) )
                    outgoingflow_per_commodity[com] = 0;

                if (edge.first == com.first
                    || edge.second == com.first) {
                    outgoingflow_per_commodity[com] += std::abs(flow_value);
                }
            }
        }

        // scale the flow values to meet one unit of flow per commodity
        for ( auto& [edge, flowMap]:f_e_st) {
            for (auto& [com, flow_value] : flowMap) {
                flow_value /= outgoingflow_per_commodity[com];
            }
        }
    }


};



#endif //OBLIVIOUSROUTING_RAECKE_FRAMEWORK_H