//
// Created by Mert Biyikli on 10.02.26.
//

#ifndef OBLIVIOUSROUTING_TREE_TRANSFORM_H
#define OBLIVIOUSROUTING_TREE_TRANSFORM_H

#include "hst.h"
#include "flat_hst.h"
#include <unordered_set>
#include <queue>
#include <cassert>

// ---------------------------------------------------------------------------
// TreeIteration<T>  —  works for T = std::shared_ptr<HSTNode>  or  T = FlatHST
// path_table[v] = precomputed path 0→v from a single SPT, built once per
// oracle call in tree_mwu and passed in here — avoids per-cut Dijkstra and
// makes all cut-paths share one consistent SPT → no cycles by construction.
// ---------------------------------------------------------------------------
template<typename T>
class TreeIteration {
public:
    T                             tree;
    std::vector<double>           distance;
    double                        lambda;


    TreeIteration(T&& _tree,
                  const std::vector<double>& _distance,
                  double _lambda)
        : tree(std::move(_tree)), distance(_distance)
        , lambda(_lambda) {}

    TreeIteration(T&& _tree,
                  std::vector<double>&& _distance,
                  double _lambda,
                  std::vector<std::vector<int>>&& _path_table)
        : tree(std::move(_tree)), distance(std::move(_distance))
        , lambda(_lambda){}

    const T&                             getTree()      const { return tree; }
    const std::vector<double>&           getDistance()  const { return distance; }
    double                               getLambda()    const { return lambda; }
};

// ---------------------------------------------------------------------------
// TreeTransform  —  single class, two overloaded distributeDemands:
//   one for std::shared_ptr<HSTNode>, one for FlatHST (flat representation).
// Both use path_table[childCenter] instead of an independent getShortestPath,
// so all cuts share one SPT root → structurally cycle-free → no removeCycles.
// ---------------------------------------------------------------------------
class TreeTransform {
public:
    const IGraph& graph;

    explicit TreeTransform(const IGraph& _graph) : graph(_graph) {}

    void transform(TreeIteration<std::shared_ptr<HSTNode>>& iter, LinearRoutingTable& table) {
        distributeDemands(iter, table);
    }
    void transform(TreeIteration<FlatHST>& iter, LinearRoutingTable& table) {
        distributeDemands(iter, table);
    }

private:
    // ------------------------------------------------------------------
    // Shared helper: apply flow along `path` (= 0 → childCenter from SPT)
    // for every (src ∈ childMembers, dst ∈ V\A) pair where src==0 or dst==0.
    // ------------------------------------------------------------------
    void applyFlow(const std::vector<int>&        path,
                   const std::vector<int>&        childMembers,
                   const std::unordered_set<int>& A,
                   double                         lambda,
                   LinearRoutingTable&            table) {
        for (int src : childMembers) {
            for (int dst : graph.getVertices()) {
                if (A.count(dst) || src == dst) continue;
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int e      = graph.getEdgeId(path[i],   path[i+1]);
                    int anti_e = graph.getEdgeId(path[i+1], path[i]);
                    assert(e != INVALID_EDGE_ID && anti_e != INVALID_EDGE_ID);
                    if (dst == 0) table.addFlow(anti_e, src, lambda);
                    if (src == 0) table.addFlow(e,      dst, lambda);
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Overload 1: pointer-based HSTNode tree
    // ------------------------------------------------------------------
    void distributeDemands(TreeIteration<std::shared_ptr<HSTNode>>& iter,
                           LinearRoutingTable& table) {
        const HSTNode*                       root       = iter.getTree().get();
        const double                         lambda     = iter.getLambda();
        const auto distance = iter.getDistance();

        std::queue<const HSTNode*> q;
        q.push(root);
        while (!q.empty()) {
            const HSTNode* current = q.front(); q.pop();
            if (!current) continue;
            for (const auto& child_ptr : current->getChildren()) {
                const HSTNode* child = child_ptr.get();
                assert(child != nullptr);
                if (child->getMembers().size() == current->getMembers().size()) {
                    q.push(child); continue;
                }
                int parentCenter = current->center;
                int childCenter  = child->center;
                assert(parentCenter != -1 && childCenter != -1);
                if (parentCenter == childCenter) { q.push(child); continue; }

                auto path = graph.getShortestPath(parentCenter, childCenter, distance);
                if (path.size() < 2) { q.push(child); continue; }

                const std::vector<int>& childMembers = child->getMembers();
                std::unordered_set<int> A(childMembers.begin(), childMembers.end());
                applyFlow(path, childMembers, A, lambda, table);
                q.push(child);
            }
        }
        removeCyclesLinear(table);
    }

    // ------------------------------------------------------------------
    // Overload 2: flat FlatHST tree
    // ------------------------------------------------------------------
    void distributeDemands(TreeIteration<FlatHST>& iter,
                           LinearRoutingTable& table) {
        const FlatHST&                           hst        = iter.getTree();
        const double                         lambda     = iter.getLambda();
        const auto distance = iter.getDistance();

        std::queue<int> q;
        q.push(hst.root());
        while (!q.empty()) {
            int cur = q.front(); q.pop();
            for (int child_idx : hst.children(cur)) {
                auto cur_members   = hst.memberRange(cur);
                auto child_members = hst.memberRange(child_idx);
                if (child_members.size() == cur_members.size()) { q.push(child_idx); continue; }

                int parentCenter = hst.nodes[cur].center;
                int childCenter  = hst.nodes[child_idx].center;
                assert(parentCenter != -1 && childCenter != -1);
                if (parentCenter == childCenter) { q.push(child_idx); continue; }

                auto path = graph.getShortestPath(parentCenter, childCenter, distance);
                if (path.size() < 2) { q.push(child_idx); continue; }

                std::vector<int> members_vec(child_members.begin(), child_members.end());
                std::unordered_set<int> A(members_vec.begin(), members_vec.end());
                applyFlow(path, members_vec, A, lambda, table);
                q.push(child_idx);
            }
        }
        removeCyclesLinear(table);
    }
public:
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
