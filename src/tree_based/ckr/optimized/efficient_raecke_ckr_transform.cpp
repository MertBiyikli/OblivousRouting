//
// Created by Mert Biyikli on 24.11.25.
//

#include "efficient_raecke_ckr_transform.h"
#include <optional>
#include <queue>
#include <unordered_set>

#include "../ckr_tree_decomposer.h"
#include <typeinfo>
#include "../../frt/frt_node.h"


/*
 * EfficientCKRTransform methods
 */



void EfficientRaeckeTransform::transform() {
    for (auto& iteration : iterations) {
        addTree(iteration);
    }
}


void EfficientRaeckeTransform::addTree(OracleTreeIteration& iteration) {
    distributeDemands(iteration.getTree(), iteration.getLambda(), iteration.getDistance());
    removeCycles();
}


const EfficientRoutingTable& EfficientRaeckeTransform::getRoutingTable() const {
    return table;
}

void EfficientRaeckeTransform::distributeDemands(const std::shared_ptr<ITreeNode> &node, double lambda, const std::vector<double>& distance) {
    if (!node) return;

    // print the tree as well
    //print_tree(node);
    std::queue<std::shared_ptr<ITreeNode>> q;
    q.push(node);

    while (!q.empty()) {
        const auto& current = q.front();
        q.pop();

        for (const auto& child : current->getChildren()) {
            if (child->getMembers().size() == current->getMembers().size()) {
                q.push(child);
            }else {
                // child and node have differnt members size -> there is a cut induced by the tree edge
                // take the representative of each cluster and push the flow induced by the cut
                // along the shortest path parent_representative -> child_representative
                int parentCenter = current->center;
                int childCenter = child->center;

                assert(parentCenter != -1
                    && childCenter != -1);

                if (parentCenter == childCenter) {
                    /*
                    // iterate through the parent until you find another center
                    std::unordered_set<int> childrenMembersSet(
                        child->getMembers().begin(),
                        child->getMembers().end()
                    );

                    int it = 0;
                    while (childrenMembersSet.find(parentCenter) != childrenMembersSet.end()
                        && (it < current->getMembers().size())){
                        parentCenter = current->getMembers()[it];
                        it++;
                    }
*/
                    q.push(child);
                    continue;
                }
                std::set<int> A;
                std::set<int> B;
                // TODO: this is nasty, but it might get the job done...
                if ( dynamic_pointer_cast<TreeNode>(node) != nullptr) {
                    A=std::set<int>(child->getMembers().begin(), child->getMembers().end());
                    // B must be relative to current, not whole graph
                    for (int v : current->getMembers()) {
                        if (!A.count(v)) B.insert(v);
                    }
                } else if
                 ( dynamic_pointer_cast<EfficientFRTTreeNode>(node) != nullptr) {
                    A = collectSubtreeVertices(child);
                    for (int v : g.getVertices()) {
                        if (!A.count(v)) B.insert(v);
                    }
                }else {
                    throw std::runtime_error("Unknown tree node type in EfficientRaeckeTransform::distributeDemands");
                }


                int fromVertex = parentCenter;
                int toVertex = childCenter;
                auto path = g.getShortestPath(parentCenter, childCenter, distance);
                if (path.size() < 2) continue;

                for (int src : A) {
                    for (int dst : B) {
                        if (src == dst) continue;
                        for (size_t i = 0; i + 1 < path.size(); ++i) {
                            std::pair<int, int> arc{path[i], path[i+1]};
                            int edgeId = g.getEdgeId(arc.first, arc.second);
                            int antiEdgeId = g.getEdgeId(arc.second, arc.first);
                            table.addFlow(edgeId, src, dst, lambda);
                            table.addFlow(antiEdgeId, dst, src, lambda);
                        }
                    }
                }
                q.push(child);
            }
        }
    }
}

std::set<int> EfficientRaeckeTransform::collectSubtreeVertices(const std::shared_ptr<ITreeNode> &node) {
    std::set<int> result(node->getMembers().begin(), node->getMembers().end());
    for (const auto& child : node->getChildren()) {
        auto sub = collectSubtreeVertices(child);
        result.insert(sub.begin(), sub.end());
    }
    return result;
}

void EfficientRaeckeTransform::normalizeOldSolutionBasedOnNewLambda(double lambda) {
    for (auto& vals : table.adj_vals) {
        for (auto& frac : vals) {
            frac *= (1.0 - lambda);
        }
    }
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

std::optional<std::vector<std::pair<int, int>>> EfficientRaeckeTransform::findCycleRec(
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
        double frac = table.getFlow(e_id, d.first, d.second);
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

std::vector<std::pair<int, int>> EfficientRaeckeTransform::findCycle(const std::pair<int, int>& d) {
    std::set<int> analyzed;
    for (int v : g.getVertices()) {
        if (analyzed.count(v)) continue;
        std::vector<int> stack;
        if (auto maybe = findCycleRec(v, analyzed, stack, d))
            return *maybe;
    }
    return {};
}