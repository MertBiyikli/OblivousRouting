//
// Created by Mert Biyikli on 24.11.25.
//

#ifndef OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_H
#define OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_H


#include "../../../solver/solver.h"
#include "../ckr_tree_decomposer.h"
#include "efficient_oracle_ckr.h"
#include "efficient_raecke_ckr_transform.h"
#include <chrono>
#include "../../../datastructures/graph_csr.h"

class EfficientRaeckeCKR : public ObliviousRoutingSolver
{
    public:
        Graph_csr* g = nullptr;

        EfficientCKR ckr_algo;
        EfficientRaeckeCKRTransform transform;
        EfficientRaeckeCKR() = default;
        ~EfficientRaeckeCKR() = default;

        void runSolve(const IGraph &g_) override {

            auto& graph = this->graphAs<Graph_csr>();     // CLEAN
            g = &graph;
            ckr_algo.setGraph(graph);


            ckr_algo.run();  // pass transform directly
            iteration_count = ckr_algo.getIterationCount();
            oracle_running_times = ckr_algo.getOracleRunningTimes();
            pure_oracle_running_times = ckr_algo.getPureOracleRunningTimes();
            scaleDownFlow(); // normalize after building
        }


        void storeFlow() override {
            auto& graph = this->graphAs<Graph_csr>();     // CLEAN
            transform.init(graph, ckr_algo.getIterations());
            transform.transform();


            // store the flow
            // given the demand map
            auto const& routingRaecke = transform.getRoutingTable();
            for (int e = 0; e < g->getNumEdges(); ++e) {
                auto& edge_flow_map = routingRaecke.adj_ids[e];
                auto& edge_flow_vals = routingRaecke.adj_vals[e];
                for (int idx = 0; idx < edge_flow_map.size(); ++idx) {
                    const auto& d = edge_flow_map[idx];
                    const auto& fraction = edge_flow_vals[idx];
                    if (d.first > d.second) continue; // skip trivial cases
                    f_e_st[{g->from[e], g->to[e]}][d] = fraction;
                }
            }

            scaleDownFlow();

        } // handled during run()

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

#endif //OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_H