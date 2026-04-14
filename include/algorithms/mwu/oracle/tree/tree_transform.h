//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_TREE_TRANSFORM_H
#define OBLIVIOUSROUTING_TREE_TRANSFORM_H


#include "../../../../data_structures/hst/pointer_hst.h"
#include "../../../../data_structures/hst/flat_hst.h"
#include <unordered_set>
#include <queue>
#include <cassert>
#include <set>

// ---------------------------------------------------------------------------
// TreeIteration<T>  —  works for T = std::shared_ptr<HSTNode>  or  T = FlatHST
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
// Both use path_cache[] instead of an independent getShortestPath,
// so all cuts share one SPT root → structurally cycle-free → no removeCycles.
// ---------------------------------------------------------------------------
class TreeTransform {
public:
    const IGraph& graph;

    explicit TreeTransform(const IGraph& _graph) : graph(_graph) {}

    void transform(TreeIteration<std::shared_ptr<HSTNode>>& iter, LinearRoutingTable& table, std::map<std::pair<int,int>, CachedPath>& path_cache) {
        buildLinearRoutingFromSpanningTree(iter, table, path_cache);
    }
    void transform(TreeIteration<FlatHST>& iter, LinearRoutingTable& table, std::map<std::pair<int,int>, CachedPath>& path_cache) {
        buildLinearRoutingFromSpanningTree(iter, table, path_cache);
    }

private:

struct UndirectedEdge {
        int u = -1;
        int v = -1;

        UndirectedEdge() = default;
        UndirectedEdge(int a, int b) {
            if (a < b) { u = a; v = b; }
            else       { u = b; v = a; }
        }

        bool operator==(const UndirectedEdge& other) const {
            return u == other.u && v == other.v;
        }
    };

    struct UndirectedEdgeHash {
        std::size_t operator()(const UndirectedEdge& e) const {
            return (static_cast<std::size_t>(e.u) << 32) ^ static_cast<std::size_t>(e.v);
        }
    };

    // ------------------------------------------------------------
    // Add all undirected edges of a path to the candidate set.
    // ------------------------------------------------------------
    void addPathEdgesToCandidateSet(
        const CachedPath& path,
        std::unordered_set<UndirectedEdge, UndirectedEdgeHash>& candidate_edges
    ) const {
        if (path.nodes.size() < 2) return;

        for (size_t i = 0; i + 1 < path.nodes.size(); ++i) {
            const int a = path.nodes[i];
            const int b = path.nodes[i + 1];
            assert(graph.getEdgeId(a, b) != INVALID_EDGE_ID);
            candidate_edges.emplace(a, b);
        }
    }

    // ------------------------------------------------------------
    // Extract a spanning tree from the candidate edge union.
    // Kruskal over the realized edges is enough here.
    // ------------------------------------------------------------
    std::vector<UndirectedEdge> extractSpanningTree(
        const std::unordered_set<UndirectedEdge, UndirectedEdgeHash>& candidate_edges
    ) const {
        std::vector<UndirectedEdge> edges(candidate_edges.begin(), candidate_edges.end());

        std::sort(edges.begin(), edges.end(), [&](const UndirectedEdge& e1, const UndirectedEdge& e2) {
            const double d1 = graph.getEdgeDistance(e1.u, e1.v);
            const double d2 = graph.getEdgeDistance(e2.u, e2.v);
            if (d1 != d2) return d1 < d2;
            if (e1.u != e2.u) return e1.u < e2.u;
            return e1.v < e2.v;
        });

        DSU dsu(graph.getNumNodes());
        std::vector<UndirectedEdge> tree_edges;
        tree_edges.reserve(graph.getNumNodes() - 1);

        for (const auto& e : edges) {
            if (dsu.unite(e.u, e.v)) {
                tree_edges.push_back(e);
                if (static_cast<int>(tree_edges.size()) == graph.getNumNodes() - 1) {
                    break;
                }
            }
        }

        return tree_edges;
    }

    // ------------------------------------------------------------
    // Root the spanning tree at root=0 and compute parent[].
    // ------------------------------------------------------------
    std::vector<int> computeParentArray(
        const std::vector<UndirectedEdge>& tree_edges,
        int root = 0
    ) const {
        const int n = graph.getNumNodes();
        std::vector<std::vector<int>> adj(n);

        for (const auto& e : tree_edges) {
            adj[e.u].push_back(e.v);
            adj[e.v].push_back(e.u);
        }

        std::vector<int> parent(n, -1);
        std::queue<int> q;
        q.push(root);
        parent[root] = root;

        while (!q.empty()) {
            const int u = q.front();
            q.pop();

            for (int v : adj[u]) {
                if (parent[v] != -1) continue;
                parent[v] = u;
                q.push(v);
            }
        }

        return parent;
    }

    // ------------------------------------------------------------
    // Convert rooted tree into the linear routing table:
    // add lambda on every edge of the unique path s -> root.
    // ------------------------------------------------------------
    void addTreePathsToLinearTable(
        const std::vector<int>& parent,
        double lambda,
        LinearRoutingTable& table,
        int root = 0
    ) const {
        const int n = graph.getNumNodes();

        for (int s = 0; s < n; ++s) {
            if (s == root) continue;
            assert(parent[s] != -1 && "Spanning tree does not reach all vertices.");

            int cur = s;
            while (cur != root) {
                const int par = parent[cur];
                assert(par != -1);

                const int e = graph.getEdgeId(cur, par);
                assert(e != INVALID_EDGE_ID);

                table.addFlow(e, s, lambda);
                cur = par;
            }
        }
    }

    // ------------------------------------------------------------
    // Pointer-HST: collect all realized parent-child center paths.
    // ------------------------------------------------------------
    void collectCandidateEdges(
        const TreeIteration<std::shared_ptr<HSTNode>>& iter,
        std::unordered_set<UndirectedEdge, UndirectedEdgeHash>& candidate_edges,
        std::map<std::pair<int,int>, CachedPath>& path_cache
    ) const {
        const HSTNode* root = iter.getTree().get();
        const auto& distance = iter.getDistance();

        std::queue<const HSTNode*> q;
        q.push(root);

        while (!q.empty()) {
            const HSTNode* current = q.front();
            q.pop();

            if (!current) continue;

            for (const auto& child_ptr : current->getChildren()) {
                const HSTNode* child = child_ptr.get();
                assert(child != nullptr);

                q.push(child);

                if (child->getMembers().size() == current->getMembers().size()) {
                    continue;
                }

                const int parentCenter = current->center;
                const int childCenter  = child->center;

                assert(parentCenter != -1 && childCenter != -1);

                if (parentCenter == childCenter) {
                    continue;
                }

                auto& path = path_cache[{parentCenter, childCenter}];
                if (path.nodes.empty()) {
                    path.nodes = graph.getShortestPathBidirectionalSearch(parentCenter, childCenter, distance);
                    // Pre-compute edge IDs for all edges on the path
                    path.edge_ids.clear();
                    for (size_t i = 0; i + 1 < path.nodes.size(); ++i) {
                        int e = graph.getEdgeId(path.nodes[i], path.nodes[i+1]);
                        path.edge_ids.push_back(e);
                    }
                }

                addPathEdgesToCandidateSet(path, candidate_edges);
            }
        }
    }

    // ------------------------------------------------------------
    // FlatHST: collect all realized parent-child center paths.
    // ------------------------------------------------------------
    void collectCandidateEdges(
        const TreeIteration<FlatHST>& iter,
        std::unordered_set<UndirectedEdge, UndirectedEdgeHash>& candidate_edges,
        std::map<std::pair<int,int>, CachedPath>& path_cache
    ) const {
        const FlatHST& hst = iter.getTree();
        const auto& distance = iter.getDistance();

        std::queue<int> q;
        q.push(hst.root());

        while (!q.empty()) {
            const int cur = q.front();
            q.pop();

            for (const int child_idx : hst.children(cur)) {
                q.push(child_idx);

                const auto cur_members   = hst.memberRange(cur);
                const auto child_members = hst.memberRange(child_idx);

                if (child_members.size() == cur_members.size()) {
                    continue;
                }

                const int parentCenter = hst.nodes[cur].center;
                const int childCenter  = hst.nodes[child_idx].center;

                assert(parentCenter != -1 && childCenter != -1);

                if (parentCenter == childCenter) {
                    continue;
                }

                auto& path = path_cache[{parentCenter, childCenter}];
                if (path.nodes.empty()) {
                    path.nodes = graph.getShortestPathBidirectionalSearch(parentCenter, childCenter, distance);
                    // Pre-compute edge IDs for all edges on the path
                    path.edge_ids.clear();
                    for (size_t i = 0; i + 1 < path.nodes.size(); ++i) {
                        int e = graph.getEdgeId(path.nodes[i], path.nodes[i+1]);
                        path.edge_ids.push_back(e);
                    }
                }
                addPathEdgesToCandidateSet(path, candidate_edges);
            }
        }
    }

    // ------------------------------------------------------------
    // Generic pipeline:
    // HST -> union of realized paths -> spanning tree -> root paths.
    // ------------------------------------------------------------
    template<typename HSTType>
    void buildLinearRoutingFromSpanningTree(
        const TreeIteration<HSTType>& iter,
        LinearRoutingTable& table,
        std::map<std::pair<int, int>, CachedPath>& path_cache
    ) const {
        std::unordered_set<UndirectedEdge, UndirectedEdgeHash> candidate_edges;
        collectCandidateEdges(iter, candidate_edges, path_cache);

        auto tree_edges = extractSpanningTree(candidate_edges);

        // Safety: in a valid connected realization we should get n-1 edges.
        // If not, something in the HST realization is disconnected.
        assert(static_cast<int>(tree_edges.size()) == graph.getNumNodes() - 1 &&
               "HST realization did not produce a connected spanning structure.");

        auto parent = computeParentArray(tree_edges, /*root=*/0);
        addTreePathsToLinearTable(parent, iter.getLambda(), table, /*root=*/0);
    }

};

#endif //OBLIVIOUSROUTING_TREE_TRANSFORM_H
