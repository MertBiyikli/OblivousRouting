//
// Created by Mert Biyikli on 10.02.26.
//

#ifndef OBLIVIOUSROUTING_TREE_TRANSFORM_H
#define OBLIVIOUSROUTING_TREE_TRANSFORM_H

#include "hst.h"
#include <set>

class TreeIteration {
public:
    std::shared_ptr<HSTNode> tree;
    std::vector<double> distance;
    double lambda;

    TreeIteration(std::shared_ptr<HSTNode> _tree, std::vector<double> _distance, double _lambda)
        : tree(std::move(_tree)), distance(std::move(_distance)), lambda(_lambda) {}

    // getter methods
    std::shared_ptr<HSTNode> getTree() const { assert(tree != nullptr); return tree; }
    std::vector<double> getDistance() const { return distance; }
    double getLambda() const { return lambda; }
};


class TreeTransform {
public:
    const IGraph& graph;

    TreeTransform(const IGraph& _graph)
        : graph(_graph){

    }


    void transform(TreeIteration& iteration, LinearRoutingTable& table) {
        distributeDemands(iteration, table);
    }

    void distributeDemands(TreeIteration& iter, LinearRoutingTable& table) {
        const auto& root = iter.getTree();
        const auto& distance = iter.getDistance();
        const auto& lambda = iter.getLambda();

        std::queue<std::shared_ptr<HSTNode>> q;
        q.push(root);

        while (!q.empty()) {
            const auto& current = q.front();
            q.pop();

            for (const auto& child : current->getChildren()) {
                if (child->getMembers().size() == current->getMembers().size()) {
                    q.push(child);
                    continue;
                }


                int parentCenter = current->center;
                int childCenter = child->center;

                assert(parentCenter != -1 && childCenter != -1);

                if (parentCenter == childCenter) {
                    q.push(child);
                    continue;
                }

                std::set<int> A;
                std::set<int> B;
                A = std::set(child->getMembers().begin(), child->getMembers().end());
                for (int v : graph.getVertices()) {
                    if (!A.count(v)) B.insert(v);
                }



                auto path = graph.getShortestPath(parentCenter, childCenter, distance);
                if (path.size() < 2) continue;

                for (int src : A) {
                    for (int dst : B) {

                        if (src == dst) continue;

                        for (size_t i = 0; i + 1 < path.size(); ++i) {
                            int e = graph.getEdgeId(path[i], path[i+1]);
                            int anti_e = graph.getEdgeId(path[i+1], path[i]);

                            assert(e != INVALID_EDGE_ID
                            && anti_e != INVALID_EDGE_ID);

                            if (dst == 0) {
                                table.addFlow(anti_e, src, lambda);
                            }
                            if (src == 0) {
                                table.addFlow(e, dst, lambda);
                            }
                        }
                    }
                }
                q.push(child);
            }
        }
        removeCyclesLinear(table);
    }

    void removeCyclesLinear(LinearRoutingTable& linearTable) {
        // Collect all sources s that have any flow in the linear table
        std::set<int> all_sources;
        for (int e = 0; e < (int)linearTable.src_ids.size(); ++e) {
            for (int s : linearTable.src_ids[e]) {
                all_sources.insert(s);
            }
        }

        for (int s : all_sources) {
            while (true) {
                auto cycle = findCycleLinear(s, linearTable);
                if (cycle.empty()) break;

                // find the minimum flow along the cycle for source s
                double minF = std::numeric_limits<double>::infinity();
                for (auto& [u, v] : cycle) {
                    int e_id = graph.getEdgeId(u, v);
                    minF = std::min(minF, linearTable.getFlow(e_id, s));
                }

                // subtract minF from every edge in the cycle
                for (auto& [u, v] : cycle) {
                    int e_id = graph.getEdgeId(u, v);
                    auto& ids  = linearTable.src_ids[e_id];
                    auto& vals = linearTable.src_flows[e_id];

                    // binary search for s in the sorted src_ids[e_id]
                    size_t len = ids.size();
                    size_t lo = 0, hi = len;
                    while (lo < hi) {
                        size_t mid = (lo + hi) >> 1;
                        if (ids[mid] < s) lo = mid + 1;
                        else              hi = mid;
                    }

                    if (lo < len && ids[lo] == s) {
                        vals[lo] -= minF;
                        if (vals[lo] <= 0.0) {
                            ids.erase( ids.begin()  + static_cast<long>(lo));
                            vals.erase(vals.begin() + static_cast<long>(lo));
                        }
                    }
                }
            }
        }
    }

    std::optional<std::vector<std::pair<int,int>>> findCycleLinearRec(
            int u,
            std::set<int>& onStack,
            std::vector<int>& stack,
            int s,
            LinearRoutingTable& linearTable) {
        if (onStack.count(u)) {
            // extract the cycle from where u first appears on the stack
            auto it = std::find(stack.begin(), stack.end(), u);
            std::vector<std::pair<int,int>> cycle;
            for (; it + 1 != stack.end(); ++it)
                cycle.emplace_back(*it, *(it + 1));
            cycle.emplace_back(stack.back(), u);
            return cycle;
        }

        stack.push_back(u);
        onStack.insert(u);

        for (int w : graph.neighbors(u)) {
            int e = graph.getEdgeId(u, w);
            if (linearTable.getFlow(e, s) <= 0.0) continue;
            if (auto res = findCycleLinearRec(w, onStack, stack, s, linearTable)) {
                return res;
            }
        }

            stack.pop_back();
            onStack.erase(u);
            return std::nullopt;
    }

    std::vector<std::pair<int,int>> findCycleLinear(int s, LinearRoutingTable& linearTable) {
        std::vector<int> stack;
        std::set<int> onStack;

        for (int v : graph.getVertices()) {
            if (onStack.count(v)) continue;
            if (auto res = findCycleLinearRec(v, onStack, stack, s, linearTable)) {
                return *res;
            }
        }
        return {};
    }
};

#endif //OBLIVIOUSROUTING_TREE_TRANSFORM_H