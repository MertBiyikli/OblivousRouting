//
// Created by Mert Biyikli on 24.11.25.
//

#ifndef OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_H
#define OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_H


#include "../../../solver/solver.h"
#include "../ckr_tree_decomposer.h"
#include "../raecke_ckr_transform.h"
#include "efficient_oracle_ckr.h"
#include <chrono>
#include "../../../datastructures/graph_csr.h"

class EfficientRaeckeCKR : public ObliviousRoutingSolver
    {


    public:
        EfficientCKR ckr_algo;
        RaeckeCKRTransform transform;
        EfficientRaeckeCKR() = default;
        ~EfficientRaeckeCKR() = default;
    void solve(const Graph &graph) override {
        return;
    }
/*
        void solve(const Graph &g) override {
            ckr_algo.debug = debug;
            ckr_algo.setGraph(g);

            ckr_algo.init(g);
            ckr_algo.run();  // pass transform directly
            iteration_count = ckr_algo.getIterationCount();
            oracle_running_times = ckr_algo.oracle_running_times;
            scaleDownFlow(); // normalize after building
        }*/

    void solve_(const Graph_csr &g) {
            ckr_algo.debug = debug;
            ckr_algo.setGraph(g);

            ckr_algo.init(g);
            ckr_algo.run();  // pass transform directly
            iteration_count = ckr_algo.getIterationCount();
            oracle_running_times = ckr_algo.oracle_running_times;
            pure_oracle_running_times = ckr_algo.pure_oracle_running_times;
            scaleDownFlow(); // normalize after building
        }

        void storeFlow() override {
            // directly add flow contribution
            for (int i = 0; i < ckr_algo.m_graphs.size(); ++i) {
                transform.addTree(ckr_algo.m_trees[i],ckr_algo.m_lambdas[i], ckr_algo.m_graphs[i]);
            }

            // store the flow
            // given the demand map
            auto const& routingRaecke = transform.getRoutingTable();
            for (const auto& [edge, demandMap] : routingRaecke) {
                for (const auto& [d, fraction] : demandMap) {
                    if (d.first > d.second) continue; // skip trivial cases
                    f_e_st[edge][d]=fraction;
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