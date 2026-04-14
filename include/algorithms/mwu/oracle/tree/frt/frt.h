//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_FRT_H
#define OBLIVIOUSROUTING_FRT_H

#include "../tree_oracle.h"

template<typename T>
class FRT : public TreeOracle<T> {
public:
    explicit FRT(IGraph& g) : TreeOracle<T>(g) {}
    explicit FRT(IGraph& g, bool mendelscaling) : TreeOracle<T>(g, mendelscaling) {}

    void computeLevelPartition(IGraph& g, HSTLevel& level, const std::vector<int>& x_perm, double delta) override {
        std::map<std::pair<int, int>, CachedPath> path_cache; // local cache for this level's shortest path computations

        level.owner.resize(g.getNumNodes());
        for (const auto& v : g.getVertices()) {
            level.owner[v] = v; // initialize owner to itself
        }
        level.R = delta;

        // compute for each node in the permutation the ball of radius delta around it and assign to the cluster of that node
        std::unordered_set<int> assigned;
        for (const auto& v : x_perm) {
            if (assigned.find(v) != assigned.end()) {
                continue;
            }

            // v is the new center for the next cluster
            level.centers.push_back(v);
            assigned.insert(v);

            for (const auto& u : g.getVertices()) {
                if (assigned.find(u) != assigned.end()) {
                    continue;
                }

                // compute shortest path distance
                auto& path = path_cache[std::make_pair(v, u)];
                if (path.nodes.empty()) {
                    // check reverse direction
                    auto& rev_path = path_cache[std::make_pair(u, v)];
                    if (!rev_path.nodes.empty()) {
                        // reverse the path and edge_ids
                        path.nodes = std::vector<int>(rev_path.nodes.rbegin(), rev_path.nodes.rend());
                        // reverse edge directions
                        std::reverse(path.nodes.begin(), path.nodes.end());

                        for (int e : path.edge_ids) {
                            path.edge_ids.push_back(g.getAntiEdge(e));
                        }
                    }else {
                        auto path_nodes = g.getShortestPathBidirectionalSearch(v, u);
                        for (size_t i = 0; i + 1 < path_nodes.size(); ++i) {
                            int from = path_nodes[i];
                            int to = path_nodes[i + 1];
                            int edge_id = g.getEdgeId(from, to);
                            if (edge_id == INVALID_EDGE_ID) {
                                throw std::runtime_error("Invalid edge in shortest path: " + std::to_string(from) + " -> " + std::to_string(to));
                            }
                            path.nodes.push_back(from);
                            path.edge_ids.push_back(edge_id);
                        }
                        path.nodes.push_back(u); // add the target node at the end
                    }
                }

                double dist = 0;
                for (size_t i = 0; i  < path.edge_ids.size(); ++i) {
                    dist += g.getEdgeDistance(path.edge_ids[i]);
                }

                if (dist <= delta) {
                    // assign u to the ball of v
                    level.owner[u] = v;
                    assigned.insert(u);
                }
            }
        }
        assert(level.centers.size() <= level.owner.size());
    }
};


#endif //OBLIVIOUSROUTING_FRT_H