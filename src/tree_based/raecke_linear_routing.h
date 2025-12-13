//
// Created by Mert Biyikli on 10.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_LINEAR_ROUTING_H
#define OBLIVIOUSROUTING_RAECKE_LINEAR_ROUTING_H

#include "../datastructures/IGraph.h"
#include "../solver/routing_table.h"
#include "raecke_oracle_iteration.h"

#include <vector>
#include <set>

class RaeckeLinearRouting {
public:
    RaeckeLinearRouting(const IGraph& _g, int _root, std::vector<OracleTreeIteration>& _iterations)
        : graph(_g), root_node(_root), iterations(_iterations) {}


    void build(LinearRoutingTable& table);

private:
    const IGraph& graph;
    int root_node;
    std::vector<OracleTreeIteration>& iterations;

    std::vector<bool> collectSubtreeVertices(const std::shared_ptr<ITreeNode>& node);
    std::set<int> _collectSubtreeVertices(const std::shared_ptr<ITreeNode>& node);
    void normalizeLambdas();
};


/*
#include <vector>
#include <memory>
#include <unordered_map>
#include <cstdint>
#include <cassert>
#include <algorithm>



class EfficientRaeckeTransformFast {
public:
    explicit EfficientRaeckeTransformFast(const IGraph& g)
        : graph(g) {}

    // Main entry:
    // Given a tree root, a demand vector b (size n), and a scalar tree_weight,
    // adds tree_weight * induced_flow(b) to `routing`.
    void addTreeAndRoute(const std::shared_ptr<ITreeNode>& root,
                         const std::vector<double>& b,
                         double tree_weight,
                         LinearRoutingTable& routing)
    {
        assert(root && "Tree root is null");
        const int n = graph.getNumNodes();
        assert((int)b.size() == n && "Demand vector has wrong size");

        FlatTree ft = buildFlatTree(root, n);

        // Bottom-up sums.
        // sum[node] = sum of b over leaves under node
        std::vector<double> sum(ft.T, 0.0);
        sum.reserve(ft.T);

        // Initialize leaves
        for (int t = 0; t < ft.T; ++t) {
            int v = ft.leaf_vertex[t];
            if (v >= 0) sum[t] = b[v];
        }

        // Postorder accumulation
        for (int x : ft.postorder) {
            int L = ft.left[x], R = ft.right[x];
            if (L >= 0) sum[x] = sum[L] + sum[R];
        }

        // Route internal-node deltas
        // delta(x) = sum(left) - sum(right)
        // route |delta| from rep(left) to rep(right), direction by sign
        for (int x = 0; x < ft.T; ++x) {
            int L = ft.left[x], R = ft.right[x];
            if (L < 0 || R < 0) continue; // leaf

            double delta = sum[L] - sum[R];
            if (delta == 0.0) continue;

            int s_rep = ft.rep[L];
            int t_rep = ft.rep[R];
            assert(s_rep >= 0 && t_rep >= 0);

            int s = (delta > 0.0 ? s_rep : t_rep);
            int t = (delta > 0.0 ? t_rep : s_rep);
            double amt = tree_weight * std::abs(delta);

            if (amt == 0.0) continue;
            routeOnGraph(s, t, amt, routing);
        }
    }

    // Optional: clear per-transform caches
    void clearCache() { path_cache.clear(); }

private:
    const IGraph& graph;

    // --- Flat tree representation for fast traversal ---
    struct FlatTree {
        int N = 0;        // #graph vertices
        int T = 0;        // #tree nodes
        int root = -1;

        std::vector<int> left, right, parent;
        std::vector<int> leaf_vertex;  // leaf -> original vertex, internal -> -1
        std::vector<int> rep;          // representative original vertex for subtree
        std::vector<int> postorder;    // internal nodes appear after children
    };

    // Map pointer identity to integer index while flattening
    struct PtrHash {
        std::size_t operator()(std::shared_ptr<ITreeNode> p) const noexcept {
            return std::hash<std::uintptr_t>{}(reinterpret_cast<std::uintptr_t>(p.get()));
        }
    };
    struct PtrEq {
        bool operator()(std::shared_ptr<ITreeNode> a, std::shared_ptr<ITreeNode> b) const noexcept { return a == b; }
    };

    FlatTree buildFlatTree(const std::shared_ptr<ITreeNode>& root, int n) {
        FlatTree ft;
        ft.N = n;

        // First pass: DFS assign ids
        std::unordered_map<std::shared_ptr<ITreeNode>, int, PtrHash, PtrEq> id;
        id.reserve(2 * n);

        std::vector<std::shared_ptr<ITreeNode>> stack;
        stack.reserve(2 * n);

        stack.push_back(root);
        id[root] = 0;

        std::vector<std::shared_ptr<ITreeNode>> order; // visitation order
        order.reserve(2 * n);

        while (!stack.empty()) {
            auto cur = stack.back();
            stack.pop_back();
            order.push_back(cur);

            auto children = cur->getChildren();
            for (auto& ch : children) {
                auto p = ch;
                if (!p) continue;
                if (!id.count(p)) {
                    int nid = (int)id.size();
                    id[p] = nid;
                    stack.push_back(p);
                }
            }
        }

        ft.T = (int)id.size();
        ft.left.assign(ft.T, -1);
        ft.right.assign(ft.T, -1);
        ft.parent.assign(ft.T, -1);
        ft.leaf_vertex.assign(ft.T, -1);
        ft.rep.assign(ft.T, -1);
        ft.postorder.clear();
        ft.postorder.reserve(ft.T);

        ft.root = id[root];

        // Second pass: fill children + detect leaves
        // We also compute postorder via explicit stack (no recursion).
        //
        // postorder stack trick: push (node, state) where state 0=enter, 1=exit
        struct Frame { std::shared_ptr<ITreeNode> node; int state; };
        std::vector<Frame> st;
        st.reserve(2 * ft.T);
        st.push_back({root, 0});
        ft.parent[ft.root] = ft.root;

        while (!st.empty()) {
            auto [node, state] = st.back();
            st.pop_back();
            int x = id[node];

            if (state == 0) {
                st.push_back({node, 1}); // exit later

                auto children = node->getChildren();
                if (children.empty()) {
                    // leaf: must be singleton member
                    const auto& mem = node->getMembers();
                    assert(mem.size() == 1 && "Leaf must have exactly one member");
                    ft.leaf_vertex[x] = mem[0];
                    ft.rep[x] = mem[0];
                } else {
                    // assume binary tree (Räcke tree). If not, you can adapt.
                    assert(children.size() == 2 && "Expected binary decomposition tree");
                    std::shared_ptr<ITreeNode> c0 = children[0];
                    std::shared_ptr<ITreeNode> c1 = children[1];
                    assert(c0 && c1);

                    int L = id[c0];
                    int R = id[c1];
                    ft.left[x] = L;
                    ft.right[x] = R;
                    ft.parent[L] = x;
                    ft.parent[R] = x;

                    st.push_back({c0, 0});
                    st.push_back({c1, 0});
                }
            } else {
                // exit: children already processed -> can compute representative
                int L = ft.left[x], R = ft.right[x];
                if (L >= 0 && R >= 0) {
                    ft.rep[x] = ft.rep[L]; // any leaf rep works; leftmost is fine
                }
                ft.postorder.push_back(x);
            }
        }

        // In postorder we included leaves too (fine). Internal nodes also included.
        return ft;
    }

    // ----------------- Shortest path + flow push -----------------

    // Optional path cache to avoid repeating s->t shortest path work within one transform run.
    // Key is ordered pair (s,t) packed into 64-bit.
    std::unordered_map<std::uint64_t, std::vector<int>> path_cache;

    static std::uint64_t packKey(int s, int t) {
        return (std::uint64_t)(std::uint32_t)s << 32 | (std::uint32_t)t;
    }

    // INTEGRATION POINT #1:
    // Return the shortest path from s to t as a list of EDGE IDS in traversal order.
    //
    // You already have an optimized getShortestPath for Graph_csr + custom PQ.
    // Hook that here (ideally returning edge ids, not vertices).
    std::vector<int> getShortestPathEdgeIds(int s, int t) {
        // TODO: replace this stub with your optimized implementation
        // Example expected behavior:
        //   return graph.getShortestPathEdgeIds(s, t);
        auto path = graph.getShortestPath(s, t);
        if (path.size() < 2) return {};
        std::vector<int> edge_ids;
        edge_ids.reserve(path.size() - 1);
        for (size_t i = 0; i + 1 < path.size(); ++i) {
            int u = path[i];
            int v = path[i + 1];
            int eId = graph.getEdgeId(u, v);
            if (eId >= 0) {
                edge_ids.push_back(eId);
            }
        }
        return edge_ids;
    }

    void routeOnGraph(int s, int t, double amt, LinearRoutingTable& routing) {
        if (s == t) return;

        const std::uint64_t key = packKey(s, t);
        auto it = path_cache.find(key);
        if (it == path_cache.end()) {
            std::vector<int> edges = getShortestPathEdgeIds(s, t);
            it = path_cache.emplace(key, std::move(edges)).first;
        }

        const std::vector<int>& path_edges = it->second;
        // If there is no path (disconnected graph), skip (or assert)
        if (path_edges.empty()) return;

        for (int eId : path_edges) {
            routing.addFlow(eId,s, amt);
        }
    }
};

*/
#endif //OBLIVIOUSROUTING_RAECKE_LINEAR_ROUTING_H