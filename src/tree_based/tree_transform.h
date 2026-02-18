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
    }
};

#endif //OBLIVIOUSROUTING_TREE_TRANSFORM_H