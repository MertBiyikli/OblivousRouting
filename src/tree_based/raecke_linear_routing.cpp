//
// Created by Mert Biyikli on 10.12.25.
//

#include "raecke_linear_routing.h"

void RaeckeLinearRouting::normalizeLambdas() {
    double lambdaSum = 0.0;
    for (const auto& iter : iterations) {
        lambdaSum += iter.getLambda();
    }

    for (auto& iter : iterations) {
        double old_lambda = iter.getLambda();
        double new_lambda = old_lambda / lambdaSum;
        iter.setLambda(new_lambda);
    }
}


std::vector<bool> RaeckeLinearRouting::collectSubtreeVertices(const std::shared_ptr<ITreeNode> &node) {
    std::vector<bool> result(graph.getNumNodes(), false);
    for (int v : node->getMembers()) {
        result[v] = true;
    }
    std::vector<bool> tmp(graph.getNumNodes(), false);
    for (const auto& child : node->getChildren()) {
        tmp = collectSubtreeVertices(child);
        for (int v = 0; v < graph.getNumNodes(); v++) {
            if (tmp[v]) result[v] = true;
        }
    }
    return result;
}


void RaeckeLinearRouting::build(LinearRoutingTable &table) {
    const int n = graph.getNumNodes();
    std::vector<int> path; // sometimes to distinct yourself from the machine, you have comment some brainrot in the code... :). FUCK I LOVE CODING!


    for (const auto& iter:iterations) {


        const double lambda = iter.getLambda();
        const auto& dist = iter.getDistance();
        auto tree = iter.getTree();

        if (!tree) continue;

        std::queue<std::shared_ptr<ITreeNode> > q;
        q.push(tree);

        while (!q.empty()) {
            auto node = q.front();
            q.pop();

            if (!node) continue;

            std::vector<bool> A;
            std::vector<bool> B;
            for (const auto& child : node->getChildren()) {
                if (child->getMembers().size() == node->getMembers().size()) {
                    q.push(child);
                    continue;
                }
                std::set<int> A = _collectSubtreeVertices(child);
                std::set<int> B;
                for (int v = 0; v<graph.getNumNodes(); v++) {
                    if (!A.count(v)) B.insert(v);
                }
                // optimized version
                /*
                                A = collectSubtreeVertices(child);
                                B.resize( graph.getNumNodes(), false);

                                for (int v = 0; v<graph.getNumNodes(); v++) {
                                    if (!A[v]) B[v]=true;
                                }

                                for (int s = 0; s<graph.getNumNodes(); s++) {
                                    if (A[s]) {
                                        for (int t = 0; t<graph.getNumNodes(); t++) {
                                            if (B[t]) {
                                                if (s == t) continue;

                                                if (s == root_node || t == root_node) {
                                                    // skip flows involving the root node
                                                    continue;
                                                }


                                                if ( true ) {
                                                    path = graph.getShortestPath(s, t, dist);
                                                    for (int i = 0; i + 1 < path.size(); ++i) {
                                                        int u = path[i];
                                                        int v = path[i + 1];

                                                        int edgeId = graph.getEdgeId(u, v);
                                                        table.addFlow(edgeId, s, lambda);
                                                        int antiEdgeId = graph.getEdgeId(v, u);
                                                        table.addFlow(antiEdgeId, t, lambda);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                */
                for (int src : A) {
                    for (int dst : B) {
                        if (src == dst) continue;


                        if (src != root_node
                            && dst == root_node) {
                            // skip flows involving the root node


                            path = graph.getShortestPath(src, dst, dist);

                            for (size_t i = 0; i + 1 < path.size(); ++i) {
                                std::pair<int, int> arc{path[i], path[i+1]};

                                int edgeId = graph.getEdgeId(arc.first, arc.second);
                                int antiEdgeId = graph.getEdgeId(arc.second, arc.first);
                                table.addFlow(edgeId, src, lambda);
                                //table.addFlow(antiEdgeId, dst, lambda);
                            }
                        }

                        if (dst != root_node
                            && src == root_node) {
                            // skip flows involving the root node


                            path = graph.getShortestPath(src, dst, dist);

                            for (size_t i = 0; i + 1 < path.size(); ++i) {
                                std::pair<int, int> arc{path[i], path[i+1]};

                                int edgeId = graph.getEdgeId(arc.first, arc.second);
                                int antiEdgeId = graph.getEdgeId(arc.second, arc.first);
                                //table.addFlow(edgeId, src, lambda);
                                table.addFlow(antiEdgeId, dst, lambda);
                            }
                        }
                    }
                }
                q.push(child);
            }
        }
    }
}

std::set<int> RaeckeLinearRouting::_collectSubtreeVertices(const std::shared_ptr<ITreeNode> &node) {
    std::set<int> result(node->getMembers().begin(), node->getMembers().end());
    for (const auto& child : node->getChildren()) {
        auto sub = collectSubtreeVertices(child);
        result.insert(sub.begin(), sub.end());
    }
    return result;
}