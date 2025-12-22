//
// Created by Mert Biyikli on 24.11.25.
//

#include "raecke_transform.h"
#include <optional>
#include <queue>
#include <unordered_set>
#include <typeinfo>
#include <__xlocale.h>


// Base class methods
void IRaeckeTransform::transform() {
    for (auto& iteration : iterations) {
        addTree(iteration);
    }
    removeCycles();
}


void IRaeckeTransform::addTree(OracleTreeIteration& iteration) {
    distributeDemands(iteration.getTree(), iteration.getLambda(), iteration.getDistance());
}

std::set<int> IRaeckeTransform::collectSubtreeVertices(const std::shared_ptr<ITreeNode> &node) {
    std::set<int> result(node->getMembers().begin(), node->getMembers().end());
    for (const auto& child : node->getChildren()) {
        auto sub = collectSubtreeVertices(child);
        result.insert(sub.begin(), sub.end());
    }
    return result;
}




const AllPairRoutingTable& EfficientRaeckeTransform::getRoutingTable() const {
    return table;
}

void EfficientRaeckeTransform::distributeDemands(const std::shared_ptr<ITreeNode> &node, double lambda, const std::vector<double>& distance) {
    if (!node) return;


    int watch_u = 2, watch_v = 3;

    // print the tree as well
    // print_tree(node);
    std::queue<std::shared_ptr<ITreeNode>> q;
    q.push(node);

    while (!q.empty()) {
        const auto& current = q.front();
        q.pop();

        for (const auto& child : current->getChildren()) {
            if (child->getMembers().size() == current->getMembers().size()) {
                q.push(child);
                continue;
            }
                // child and node have differnt members size -> there is a cut induced by the tree edge
                // take the representative of each cluster and push the flow induced by the cut
                // along the shortest path parent_representative -> child_representative
                int parentCenter = current->center;
                int childCenter = child->center;

                assert(parentCenter != -1
                    && childCenter != -1);

                if (parentCenter == childCenter) {
                    q.push(child);
                    continue;
                }
                std::set<int> A;
                std::set<int> B;
                A = std::set(child->getMembers().begin(), child->getMembers().end());
                for (int v : g.getVertices()) {
                    if (!A.count(v)) B.insert(v);
                }



                auto path = g.getShortestPath(parentCenter, childCenter, distance);
                if (path.size() < 2) continue;

                for (int src : A) {
                    for (int dst : B) {

                        if (src == dst) continue;
/*
                        if( (src == watch_u && dst == watch_v)
                        || (src == watch_v && dst == watch_u)) {
                            std::cout << "Child member: ";
                            for (auto v : A) {
                                std::cout << v << " ";
                            }
                            std::cout << "\n";
                            std::cout << "Parent member: ";
                            for (auto v : B) {
                                std::cout << v << " ";
                            }
                            std::cout << "\n";
                            std::cout << "Pushing flow for demand pair (" << src << ", " << dst << ") along path: ";
                            for (auto v : path) {
                                std::cout << v << " ";
                            }
                            std::cout << "\n";
                            std::cout << "Lambda : " << lambda << "\n";
                        }*/
                        for (size_t i = 0; i + 1 < path.size(); ++i) {

                            int e = g.getEdgeId(path[i], path[i+1]);
                            int anti_e = g.getEdgeId(path[i+1], path[i]);

                            assert(e != INVALID_EDGE_ID && anti_e != INVALID_EDGE_ID);
                            table.addFlow(anti_e, src, dst, lambda);
                            table.addFlow(e, dst, src, lambda);
/*
                            auto& forward = edge2demand2flow[e];
                            auto& reverse = edge2demand2flow[anti_e];


                            if (forward.count({src, dst}) == 0) {
                                forward[{src, dst}] = 0.0;
                            }
                            if (reverse.count({dst, src}) == 0) {
                                reverse[{dst, src}] = 0.0;
                            }

                            forward[{src, dst}] += lambda;
                            reverse[{dst, src}] += lambda;
*/
                        }
                    }
                }
                q.push(child);
        }
    }
}


bool EfficientRaeckeTransform::routingTableIsValid() {
    // check the flow conservation and net flow = 1 for each commodity
    // check the map DS
    std::map<std::pair<int, int>, double> net_flow;
    for (int e = 0; e<g.getNumEdges(); e++) {
        const auto& demand2flow = edge2demand2flow[e];
        for (const auto& [d, f] : demand2flow) {
            if (d.first == g.edgeEndpoints(e).first
                || d.first == g.edgeEndpoints(e).second) {
                net_flow[d] += f;
            }
        }
    }
    bool all_unit_flow = true;
    for (const auto& [d, f] : net_flow) {
        if (std::abs(f - 1.0) > SOFT_EPS) {
            all_unit_flow = false;
            std::cout << "Commodity " << d.first << " -> " << d.second << " has net flow "
                      << f << ")\n";
        }
    }
    return all_unit_flow;
}





void EfficientRaeckeTransform::removeCycles() {
    std::set<std::pair<int, int>> all;
    for (int e = 0; e<table.adj_ids.size(); e++) {
        const auto& ids  = table.adj_ids[e];
        for (const auto& id : ids) {
            all.insert(id);
        }
    }

    for (auto d : all) {
        while (true) {
            bool debug = false;
            if (d.first == 2 && d.second == 3) {
                debug = true;
            }
            auto cycle = findCycle(d);
            if (cycle.empty()) break;

            double minF = std::numeric_limits<double>::infinity();
            for (auto e : cycle) {
                int e_id = g.getEdgeId(e.first, e.second);
                double frac = table.getFlow(e_id, d.first, d.second);
                minF = std::min(minF,frac);
            }

            for (auto e : cycle) {
                int e_id = g.getEdgeId(e.first, e.second);
                auto& dmap = table.adj_vals[e_id];
                // find index of (d.first, d.second) in table.adj_ids[e_id]
                const auto& ids = table.adj_ids[e_id];
                size_t len = ids.size();
                size_t lo = 0, hi = len;
                while (lo < hi) {
                    const size_t mid = (lo + hi) >> 1;
                    const auto& mid_val = ids[mid];
                    if (mid_val < std::make_pair(d.first, d.second))
                        lo = mid + 1;
                    else
                        hi = mid;
                }

                // if found, decrease the fraction
                if (lo < len && ids[lo] == std::make_pair(d.first, d.second)) {
                    dmap[lo] -= minF;
                    if (dmap[lo] <= 0) {
                        // remove the entry
                        dmap.erase(dmap.begin() + static_cast<long>(lo));
                        table.adj_ids[e_id].erase(table.adj_ids[e_id].begin() + static_cast<long>(lo));
                    }
                }

                // also erase for the anti edge
                int anti_e_id = g.getEdgeId(e.second, e.first);
                auto& rmap = table.adj_vals[anti_e_id];
                const auto& rids = table.adj_ids[anti_e_id];
                size_t rlen = rids.size();
                size_t rlo = 0, rhi = rlen;
                while (rlo < rhi) {
                    const size_t mid = (rlo + rhi) >> 1;
                    const auto& mid_val = rids[mid];
                    if (mid_val < std::make_pair(d.second, d.first))
                        rlo = mid + 1;
                    else
                        rhi = mid;
                }

                // if found, decrease the fraction
                if (rlo < rlen && rids[rlo] == std::make_pair(d.second, d.first)) {
                    rmap[rlo] -= minF;
                    if (rmap[rlo] <= 0) {
                        // remove the entry
                        rmap.erase(rmap.begin() + static_cast<long>(rlo));
                        table.adj_ids[anti_e_id].erase(table.adj_ids[anti_e_id].begin() + static_cast<long>(rlo));
                    }
                }
            }

        }
    }
}

std::optional<std::vector<std::pair<int,int>>>
EfficientRaeckeTransform::findCycleRec(
    int u,
    std::set<int>& onStack,
    std::vector<int>& stack,
    const std::pair<int,int>& d
) {
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

    for (int w : g.neighbors(u)) {
        int e = g.getEdgeId(u, w);
        if (table.getFlow(e, d.first, d.second) <= 0) continue;

        if (auto res = findCycleRec(w, onStack, stack, d))
            return res;
    }

    stack.pop_back();
    onStack.erase(u);
    return std::nullopt;
}


std::vector<std::pair<int,int>>
EfficientRaeckeTransform::findCycle(const std::pair<int,int>& d) {
    std::vector<int> stack;
    std::set<int> onStack;

    for (int v : g.getVertices()) {
        if (onStack.count(v)) continue;
        if (auto res = findCycleRec(v, onStack,stack, d))
            return *res;
    }
    return {};
}






const LinearRoutingTable& LinearEfficientRaeckeTransform::getRoutingTable() const {
    return linear_tb;
}

void LinearEfficientRaeckeTransform::distributeDemands(const std::shared_ptr<ITreeNode> &tree, double lambda, const std::vector<double>& distance)  {
    if (!tree) return;

    // print the tree as well
    //print_tree(tree);
    std::queue<std::shared_ptr<ITreeNode>> q;
    for (const auto& child : tree->getChildren()) {
        q.push(child);
    }

    while (!q.empty()) {
        const auto& current = q.front();
        const auto& parent = current->getParent();
        assert(current != nullptr);
        assert(parent != nullptr);
        q.pop();




        // child and node have different members size -> there is a cut induced by the tree edge
        // take the representative of each cluster and push the flow induced by the cut
        // along the shortest path parent_representative -> child_representative
        int parentCenter = parent->center;
        int childCenter = current->center;

        assert(parentCenter != -1
            && childCenter != -1);

        if (parentCenter == childCenter) {
            for (const auto& child : current->getChildren()) {
                q.push(child);
            }
            continue;
        }
        // IMPORTANT: A must be the CHILD SUBTREE (as in your Efficient version)
        // Using only child->getMembers() is wrong if members are not subtree-closed.
        std::set<int> A = collectSubtreeVertices(current);
        if (A.empty()) continue;

        // print A
        std::cout << "Subtree vertices (A) for child center " << childCenter << ": ";
        for (const auto& v : A) {
            std::cout << v << " ";
        }
        std::cout << std::endl;

        // shortest path between the representatives in the ORIGINAL graph metric
        auto path = g.getShortestPath(childCenter, parentCenter, distance);
        if (path.size() < 2) continue;

        // print path
        std::cout << "Shortest path from " << childCenter << " to " << parentCenter << ": ";
        for (const auto& v : path) {
            std::cout << v << " ";
        }
        std::cout << std::endl;

        // Update the linear operator: for each u in A, add +lambda along the path
        for (int u : A) {
            if ( u == root) continue;
            for (size_t i = 0; i + 1 < path.size(); ++i) {
                const int a = path[i];
                const int b = path[i + 1];
                const int e = g.getEdgeId(a, b);
                if (e < 0) continue;

                std::cout << "Pushing flow along edge: " << a << " -> " << b
                          << " for demand to " << u << " with lambda " << lambda << "\n";

                // Column update: M[e, u] += lambda
                linear_tb.addFlow(e, u, lambda);
            }
        }

        for (const auto& child : current->getChildren()) {
            q.push(child);
        }
    }
    // linear_tb.printFlows(g);
    removeCycles();
}




void LinearEfficientRaeckeTransform::removeCycles() {
    std::set<std::pair<int, int>> all;
    for (int e = 0; e<linear_tb.src_ids.size(); e++) {
        const auto& ids  = linear_tb.src_ids[e];
        for (const auto& id : ids) {
            all.insert({id, root});
        }
    }

    for (auto d : all) {
        while (true) {
            auto cycle = findCycle(d);
            if (cycle.empty()) break;

            double minF = std::numeric_limits<double>::infinity();
            for (auto e : cycle) {
                int e_id = g.getEdgeId(e.first, e.second);
                double frac = linear_tb.getFlow(e_id, d.first);
                minF = std::min(minF,frac);
            }

            for (auto e : cycle) {
                int e_id = g.getEdgeId(e.first, e.second);
                auto& dmap = linear_tb.src_flows[e_id];
                // find index of (d.first, d.second) in table.adj_ids[e_id]
                const auto& ids = linear_tb.src_ids[e_id];
                size_t len = ids.size();
                size_t lo = 0, hi = len;
                while (lo < hi) {
                    const size_t mid = (lo + hi) >> 1;
                    const auto& mid_val = ids[mid];
                    if (mid_val < d.first)
                        lo = mid + 1;
                    else
                        hi = mid;
                }

                // if found, decrease the fraction
                if (lo < len && ids[lo] == d.first) {
                    dmap[lo] -= minF;
                    if (dmap[lo] <= 0) {
                        // remove the entry
                        dmap.erase(dmap.begin() + static_cast<long>(lo));
                        linear_tb.src_ids[e_id].erase(linear_tb.src_ids[e_id].begin() + static_cast<long>(lo));
                    }
                }

            }

        }
    }
}

std::optional<std::vector<std::pair<int, int>>> LinearEfficientRaeckeTransform::findCycleRec(
        int u,
        std::set<int>& analyzed,
        std::vector<int>& stack,
        const std::pair<int, int>& d) {
    if (std::find(stack.begin(), stack.end(), u) != stack.end())
        return std::vector<std::pair<int, int>>{};
    if (analyzed.count(u))
        return std::nullopt;

    analyzed.insert(u);
    stack.push_back(u);

    for (int w : g.neighbors(u)) {
        int e_id = g.getEdgeId(u, w);
        double frac = linear_tb.getFlow(e_id, d.first);
        if (frac <= 0) continue;

        if (auto child = findCycleRec(w, analyzed, stack, d)) {
            auto cycle = *child;
            cycle.insert(cycle.begin(), {u, w});
            return cycle;
        }
    }

    stack.pop_back();
    return std::nullopt;
}

std::vector<std::pair<int, int>> LinearEfficientRaeckeTransform::findCycle(const std::pair<int, int>& d) {
    std::set<int> analyzed;
    for (int v : g.getVertices()) {
        if (analyzed.count(v)) continue;
        std::vector<int> stack;
        if (auto maybe = findCycleRec(v, analyzed, stack, d))
            return *maybe;
    }
    return {};
}