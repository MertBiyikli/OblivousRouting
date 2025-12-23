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
    // removeCycles();
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

bool EfficientRaeckeTransform::routingTableIsValid() const {
    return table.isValid(g);
}

void EfficientRaeckeTransform::distributeDemands(const std::shared_ptr<ITreeNode> &node, double lambda, const std::vector<double>& distance) {
    if (!node) return;
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

                        for (size_t i = 0; i + 1 < path.size(); ++i) {

                            int e = g.getEdgeId(path[i], path[i+1]);
                            int anti_e = g.getEdgeId(path[i+1], path[i]);

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






void EfficientRaeckeTransform::removeCycles() {
    std::set<std::pair<int, int> > all;
    for (int e = 0; e<table.adj_ids.size(); e++) {
        const auto& ids  = table.adj_ids[e];
        for (const auto& id : ids) {
            int s = id / g.getNumNodes();
            int t = id % g.getNumNodes();
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
                int e_id = g.getEdgeId(e.first, e.second);
                double frac = (table.getFlow(e_id, d.first, d.second));
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
                    if (mid_val < getCommodityID(g.getNumNodes(), d.first, d.second))
                        lo = mid + 1;
                    else
                        hi = mid;
                }

                // if found, decrease the fraction
                if (lo < len && ids[lo] == getCommodityID(g.getNumNodes(), d.first, d.second)) {
                    dmap[lo] -= minF;
                    if (dmap[lo] <= 0) {
                        // remove the entry
                        dmap.erase(dmap.begin() + static_cast<long>(lo));
                        table.adj_ids[e_id].erase(table.adj_ids[e_id].begin() + static_cast<long>(lo));
                    }
                }
/*
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
                }*/
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

void LinearEfficientRaeckeTransform::distributeDemands(const std::shared_ptr<ITreeNode> &node, double lambda, const std::vector<double>& distance)  {
    if (!node) return;
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
                        // we only care about the flow from src to root
                        if (dst == root) {
                            for (size_t i = 0; i + 1 < path.size(); ++i) {

                                // int e = g.getEdgeId(path[i], path[i+1]);
                                int anti_e = g.getEdgeId(path[i+1], path[i]);

                                assert(anti_e != INVALID_EDGE_ID);
                                linear_tb.addFlow(anti_e, src, lambda);
                                //linear_tb.addFlow(e, dst, src, lambda);

                            }
                        }
/*
                        if ( src == root) {
                            for (size_t i = 0; i + 1 < path.size(); ++i) {

                                int e = g.getEdgeId(path[i], path[i+1]);
                                // int anti_e = g.getEdgeId(path[i+1], path[i]);

                                assert(e != INVALID_EDGE_ID);
                                linear_tb.addFlow(e, dst, lambda);
                                //linear_tb.addFlow(e, dst, src, lambda);

                            }
                        }*/
                    }
                }
                q.push(child);
        }
    }
}




void LinearEfficientRaeckeTransform::removeCycles() {
    std::set<int> all;
    for (int e = 0; e<linear_tb.src_ids.size(); e++) {
        const auto& ids  = linear_tb.src_ids[e];
        for (const auto& id : ids) {
            all.insert(id);
        }
    }

    for (auto d : all) {
        while (true) {
            bool debug = false;

            auto cycle = findCycle(d);
            if (cycle.empty()) break;

            double minF = std::numeric_limits<double>::infinity();
            for (auto e : cycle) {
                int e_id = g.getEdgeId(e.first, e.second);
                double frac = linear_tb.getFlow(e_id, d);
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
                    if (mid_val < d)
                        lo = mid + 1;
                    else
                        hi = mid;
                }

                // if found, decrease the fraction
                if (lo < len && ids[lo] == d) {
                    dmap[lo] -= minF;
                    if (dmap[lo] <= 0) {
                        // remove the entry
                        dmap.erase(dmap.begin() + static_cast<long>(lo));
                        linear_tb.src_ids[e_id].erase(linear_tb.src_ids[e_id].begin() + static_cast<long>(lo));
                    }
                }

                // also erase for the anti edge
                int anti_e_id = g.getEdgeId(e.second, e.first);
                auto& rmap = linear_tb.src_flows[anti_e_id];
                const auto& rids = linear_tb.src_ids[anti_e_id];
                size_t rlen = rids.size();
                size_t rlo = 0, rhi = rlen;
                while (rlo < rhi) {
                    const size_t mid = (rlo + rhi) >> 1;
                    const auto& mid_val = rids[mid];
                    if (mid_val < d)
                        rlo = mid + 1;
                    else
                        rhi = mid;
                }

                // if found, decrease the fraction
                if (rlo < rlen && rids[rlo] == d) {
                    rmap[rlo] -= minF;
                    if (rmap[rlo] <= 0) {
                        // remove the entry
                        rmap.erase(rmap.begin() + static_cast<long>(rlo));
                        linear_tb.src_ids[anti_e_id].erase(linear_tb.src_ids[anti_e_id].begin() + static_cast<long>(rlo));
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
        const int& d) {
    if (std::find(stack.begin(), stack.end(), u) != stack.end())
        return {};
    if (analyzed.count(u))
        return std::nullopt;

    analyzed.insert(u);
    stack.push_back(u);

    for (int w : g.neighbors(u)) {
        int e_id = g.getEdgeId(u, w);
        double frac = linear_tb.getFlow(e_id, d);
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

std::vector<std::pair<int, int>> LinearEfficientRaeckeTransform::findCycle(const int& d) {
    std::set<int> analyzed;
    for (int v : g.getVertices()) {
        if (analyzed.count(v)) continue;
        std::vector<int> stack;
        if (auto maybe = findCycleRec(v, analyzed, stack, d))
            return *maybe;
    }
    return {};
}