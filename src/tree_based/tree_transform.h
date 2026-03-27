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
#include <set>

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
                  double _lambda)
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

        bool flowGraphIsAcyclic = true;

        for (int s : all_sources) {
            while (true) {
                auto cycle = findCycleLinear(s, linearTable);
                if (cycle.empty()) {
                    flowGraphIsAcyclic &= true;
                }else {
                    flowGraphIsAcyclic = false;
                }
                if (cycle.empty()) break;

                std::cout << "Found cycle with source s: " << s << std::endl;
                for (auto & [u, v] : cycle) {
                    std::cout << "( " << u << " , " << v << ") ";
                }
                std::cout << std::endl;

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
        if (flowGraphIsAcyclic) {
            std::cout << "The flow graph is already acyclic. No cycles were found." << std::endl;
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

    // ------------------------------------------------------------------
// Tarjan SCC on the subgraph of edges with positive flow for source s.
// Returns a list of SCCs; each SCC with size >= 2 contains a cycle.
// Uses iterative DFS to avoid stack overflow on large graphs.
// ------------------------------------------------------------------
std::vector<std::vector<int>> tarjanSCC(int s, LinearRoutingTable& linearTable) {
    const std::vector<int>& verts = graph.getVertices();
    int N = *std::max_element(verts.begin(), verts.end()) + 1;

    std::vector<int>  index_of(N, -1);
    std::vector<int>  lowlink(N, 0);
    std::vector<bool> on_stack(N, false);
    std::vector<int>  stk;
    std::vector<std::vector<int>> sccs;
    int idx = 0;

    struct Frame { int v; int ni; };
    std::vector<Frame> call_stack;

    auto strongconnect = [&](int root) {
        call_stack.push_back({root, 0});
        index_of[root] = lowlink[root] = idx++;
        stk.push_back(root);
        on_stack[root] = true;

        while (!call_stack.empty()) {
            auto& [v, ni] = call_stack.back();
            std::vector<int> nbrs;
            for (auto u : graph.neighbors(v)) {
                nbrs.push_back(u);
            }

            bool pushed = false;
            while (ni < (int)nbrs.size()) {
                int w = nbrs[ni];
                int e = graph.getEdgeId(v, w);
                ++ni;
                if (linearTable.getFlow(e, s) <= 0.0) continue;

                if (index_of[w] == -1) {
                    // not yet visited — recurse
                    index_of[w] = lowlink[w] = idx++;
                    stk.push_back(w);
                    on_stack[w] = true;
                    call_stack.push_back({w, 0});
                    pushed = true;
                    break;
                } else if (on_stack[w]) {
                    lowlink[v] = std::min(lowlink[v], index_of[w]);
                }
            }

            if (!pushed) {
                // All neighbours processed — pop frame
                call_stack.pop_back();
                if (!call_stack.empty()) {
                    int parent = call_stack.back().v;
                    lowlink[parent] = std::min(lowlink[parent], lowlink[v]);
                }
                // Root of an SCC?
                if (lowlink[v] == index_of[v]) {
                    std::vector<int> scc;
                    while (true) {
                        int w = stk.back(); stk.pop_back();
                        on_stack[w] = false;
                        scc.push_back(w);
                        if (w == v) break;
                    }
                    sccs.push_back(std::move(scc));
                }
            }
        }
    };

    for (int v : verts) {
        if (index_of[v] == -1) strongconnect(v);
    }
    return sccs;
}

// ------------------------------------------------------------------
// Find one cycle within a given SCC using iterative DFS.
// Only traverses edges with positive flow for source s,
// restricted to vertices inside the SCC.
// ------------------------------------------------------------------
std::vector<std::pair<int,int>> findCycleInSCC(
        const std::vector<int>&  scc_verts,
        int                      s,
        LinearRoutingTable&      linearTable)
{
    std::unordered_set<int> in_scc(scc_verts.begin(), scc_verts.end());
    std::unordered_map<int,int>  parent;
    std::unordered_set<int>      visited;
    std::unordered_set<int>      on_path;
    std::vector<int>             dfs_stack;

    for (int root : scc_verts) {
        if (visited.count(root)) continue;
        dfs_stack.push_back(root);

        while (!dfs_stack.empty()) {
            int v = dfs_stack.back();
            if (!visited.count(v)) {
                visited.insert(v);
                on_path.insert(v);
            }

            bool pushed = false;
            for (int w : graph.neighbors(v)) {
                if (!in_scc.count(w)) continue;
                int e = graph.getEdgeId(v, w);
                if (linearTable.getFlow(e, s) <= 0.0) continue;

                if (on_path.count(w)) {
                    // Found back-edge v->w — reconstruct the cycle w -> ... -> v -> w
                    std::vector<std::pair<int,int>> cycle;
                    cycle.emplace_back(v, w);
                    int cur = v;
                    while (cur != w) {
                        int par = parent.at(cur);
                        cycle.emplace_back(par, cur);
                        cur = par;
                    }
                    std::reverse(cycle.begin(), cycle.end());
                    return cycle;
                }
                if (!visited.count(w)) {
                    parent[w] = v;
                    dfs_stack.push_back(w);
                    pushed = true;
                    break;
                }
            }

            if (!pushed) {
                on_path.erase(v);
                dfs_stack.pop_back();
            }
        }
    }
    return {};
}

// ------------------------------------------------------------------
// Subtract minF of a cycle for source s from linearTable.
// Uses binary search since src_ids[e] is kept sorted.
// ------------------------------------------------------------------
void subtractCycle(const std::vector<std::pair<int,int>>& cycle,
                   int s, double minF,
                   LinearRoutingTable& linearTable)
{
    for (auto& [u, v] : cycle) {
        int e_id   = graph.getEdgeId(u, v);
        auto& ids  = linearTable.src_ids[e_id];
        auto& vals = linearTable.src_flows[e_id];

        size_t lo = 0, hi = ids.size();
        while (lo < hi) {
            size_t mid = (lo + hi) >> 1;
            if (ids[mid] < s) lo = mid + 1;
            else              hi = mid;
        }
        if (lo < ids.size() && ids[lo] == s) {
            vals[lo] -= minF;
            if (vals[lo] <= 0.0) {
                ids.erase( ids.begin()  + static_cast<long>(lo));
                vals.erase(vals.begin() + static_cast<long>(lo));
            }
        }
    }
}

// ------------------------------------------------------------------
// Main entry: remove all cycles using Tarjan SCCs.
// For each source s:
//   1. Compute SCCs on the positive-flow subgraph.
//   2. For every non-trivial SCC (size >= 2):
//      a. Find a cycle within the SCC.
//      b. Subtract its minimum flow (cancels at least one edge).
//      c. Repeat until the SCC contains no more cycles.
//   3. Recompute SCCs if any edge was zeroed (changed == true).
// ------------------------------------------------------------------
void removeCyclesUsingTarjanSCCAlgorithm(LinearRoutingTable& linearTable) {
    std::set<int> all_sources;
    for (int e = 0; e < (int)linearTable.src_ids.size(); ++e)
        for (int s : linearTable.src_ids[e])
            all_sources.insert(s);

    for (int s : all_sources) {
        bool changed = true;
        while (changed) {
            changed = false;
            auto sccs = tarjanSCC(s, linearTable);
            for (auto& scc : sccs) {
                if (scc.size() < 2) continue;
                // eliminate all cycles within this SCC before moving on
                while (true) {
                    auto cycle = findCycleInSCC(scc, s, linearTable);
                    if (cycle.empty()) break;
                    double minF = std::numeric_limits<double>::infinity();
                    for (auto& [u, v] : cycle) {
                        int e_id = graph.getEdgeId(u, v);
                        minF = std::min(minF, linearTable.getFlow(e_id, s));
                    }
                    subtractCycle(cycle, s, minF, linearTable);
                    changed = true;
                }
            }
        }
    }
}

};

#endif //OBLIVIOUSROUTING_TREE_TRANSFORM_H
