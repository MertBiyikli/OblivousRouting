//
// Created by Mert Biyikli on 12.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_MWU_RANDOM_H
#define OBLIVIOUSROUTING_RAECKE_MWU_RANDOM_H

#include "../raecke_mwu.h"
#include "raecke_oracle_random.h"

class RaeckeMWU_Random : public RaeckeMWU {
public:
    RaeckeMWU_Random(IGraph& g, int root)
        : RaeckeMWU(g, root, std::make_unique<RandomOracle>(g)) {}

    virtual void transformSolution(LinearRoutingTable& table) override {
        // Implement the transformation logic specific to the Random Oracle here
         EfficientRaeckeTransform transform(graph, iteration);
        transform.transform();
        auto tree_table = transform.getRoutingTable();

        const int m = graph.getNumEdges();
        const int n = graph.getNumNodes();
        const int R = root;

        // 3) First pass: compute total outgoing flow for each s → root
        std::vector<double> out(n, 0.0);

        for (int e = 0; e < m; ++e) {
            auto [u, v] = graph.edgeEndpoints(e);
            const auto& ids  = tree_table.adj_ids[e];
            const auto& vals = tree_table.adj_vals[e];

            for (int k = 0; k < (int)ids.size(); ++k) {
                int s = ids[k].first;
                int t = ids[k].second;
                if (t != R) continue;          // only s → root commodities

                double f = vals[k];
                // We consider flow "leaving s" if edge incident to s
                if (u == s || v == s) {
                    out[s] += std::fabs(f);
                }
            }
        }

        // 4) Second pass: normalize to 1 unit and store basis flows
        for (int e = 0; e < m; ++e) {
            const auto& ids  = tree_table.adj_ids[e];
            const auto& vals = tree_table.adj_vals[e];

            for (int k = 0; k < (int)ids.size(); ++k) {
                int s = ids[k].first;
                int t = ids[k].second;
                if (t != R) continue;          // only s → root

                if (out[s] <= 1e-15) continue; // degenerate / isolated source

                double flow_sR = vals[k] / out[s];  // now 1-unit normalization
                if (std::fabs(flow_sR) < 1e-15) continue;

                // This writes f_e(s, root) into the LinearRoutingTable
                table.addFlow(e, s, flow_sR);
            }
        }

        // print out tree table
        for (int e = 0; e < graph.getNumEdges(); ++e) {
            for (int j = 0; j < tree_table.adj_ids[e].size(); ++j) {
                int source = tree_table.adj_ids[e][j].first;
                int target = tree_table.adj_ids[e][j].second;
                double flow = tree_table.adj_vals[e][j];
                if ( source == 0 && target == 1) {
                    if (flow != 0.0) {
                        auto [head, tail] = graph.edgeEndpoints(e);
                       //  std::cout << "Edge " << head << " / " << tail << " has flow " << flow << " for source " << source << " to target " << target << "\n";
                    }
                }
            }
        }
    }
};

#endif //OBLIVIOUSROUTING_RAECKE_MWU_RANDOM_H