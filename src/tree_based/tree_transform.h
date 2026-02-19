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

    // setter methods
    void setTree(std::shared_ptr<HSTNode> _tree) { tree = std::move(_tree); }
    void setDistance(std::vector<double> _distance) { distance = std::move(_distance); }
    void setLambda(double _lambda) { lambda = _lambda; }

    // getter methods
    std::shared_ptr<HSTNode> getTree() const { assert(tree != nullptr); return tree; }
    std::vector<double> getDistance() const { return distance; }
    double getLambda() const { return lambda; }
};


class TreeTransform {
public:
    const IGraph& graph;
    std::vector<TreeIteration>& iterations;

    AllPairRoutingTable table;

    TreeTransform(const IGraph& _graph, std::vector<TreeIteration>& _iterations)
        : graph(_graph), iterations(_iterations) {
        table.init(graph);
    }

    void transform() {
        for (auto& iteration : iterations) {
            distributeDemands(iteration);
        }
    }

    const AllPairRoutingTable& getRoutingTable() const {
        return table;
    }

    void distributeDemands(TreeIteration& iter) {
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
                            table.addFlow(anti_e, src, dst, lambda);
                            table.addFlow(e, dst, src, lambda);

                        }
                    }
                }
                q.push(child);
            }
        }
        removeCycles();
    }

    void removeCycles() {
        std::set<std::pair<int, int> > all;
        for (int e = 0; e<table.adj_ids.size(); e++) {
            const auto& ids  = table.adj_ids[e];
            for (const auto& id : ids) {
                int s = id / graph.getNumNodes();
                int t = id % graph.getNumNodes();
                all.insert({s, t});
            }
        }

        for (auto d : all) {
            while (true) {
                bool debug = false;

                auto cycle = findCycle(d);
                if (cycle.empty()) break;

                double minF = std::numeric_limits<double>::infinity();
                for (auto e : cycle) {
                    int e_id = graph.getEdgeId(e.first, e.second);
                    double frac = (table.getFlow(e_id, d.first, d.second));
                    minF = std::min(minF,frac);
                }

                for (auto e : cycle) {
                    int e_id = graph.getEdgeId(e.first, e.second);
                    auto& dmap = table.adj_vals[e_id];
                    // find index of (d.first, d.second) in table.adj_ids[e_id]
                    const auto& ids = table.adj_ids[e_id];
                    size_t len = ids.size();
                    size_t lo = 0, hi = len;
                    while (lo < hi) {
                        const size_t mid = (lo + hi) >> 1;
                        const auto& mid_val = ids[mid];
                        if (mid_val < getCommodityID(graph.getNumNodes(), d.first, d.second))
                            lo = mid + 1;
                        else
                            hi = mid;
                    }

                    // if found, decrease the fraction
                    if (lo < len && ids[lo] == getCommodityID(graph.getNumNodes(), d.first, d.second)) {
                        dmap[lo] -= minF;
                        if (dmap[lo] <= 0) {
                            // remove the entry
                            dmap.erase(dmap.begin() + static_cast<long>(lo));
                            table.adj_ids[e_id].erase(table.adj_ids[e_id].begin() + static_cast<long>(lo));
                        }
                    }

                }

            }
        }
    }

    std::optional<std::vector<std::pair<int, int>>> findCycleRec(
        int u,
        std::set<int>& onStack,
        std::vector<int>& stack,
        const std::pair<int, int>& d) {
        // Found a back-edge → cycle
        if (onStack.count(u)) {
            // extract the cycle vertices
            auto it = std::find(stack.begin(), stack.end(), u);
            std::vector<std::pair<int,int>> cycle;
            for (; it + 1 != stack.end(); ++it) {
                cycle.emplace_back(*it, *(it + 1));
            }
            cycle.emplace_back(stack.back(), u); // close the cycle
            return cycle;
        }

        stack.push_back(u);
        onStack.insert(u);

        for (int w : graph.neighbors(u)) {
            int e = graph.getEdgeId(u, w);
            if (table.getFlow(e, d.first, d.second) <= 0) continue;

            if (auto res = findCycleRec(w, onStack, stack, d))
                return res;
        }

        stack.pop_back();
        onStack.erase(u);
        return std::nullopt;
    }

    std::vector<std::pair<int, int>> findCycle(const std::pair<int, int>& d) {
        std::vector<int> stack;
        std::set<int> onStack;

        for (int v : graph.getVertices()) {
            if (onStack.count(v)) continue;
            if (auto res = findCycleRec(v, onStack,stack, d))
                return *res;
        }
        return {};
    }
};

#endif //OBLIVIOUSROUTING_TREE_TRANSFORM_H