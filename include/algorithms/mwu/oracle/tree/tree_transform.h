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
// Both use path_table[childCenter] instead of an independent getShortestPath,
// so all cuts share one SPT root → structurally cycle-free → no removeCycles.
// ---------------------------------------------------------------------------
class TreeTransform {
public:
    const IGraph& graph;

    explicit TreeTransform(const IGraph& _graph) : graph(_graph) {}

    void transform(TreeIteration<std::shared_ptr<HSTNode>>& iter, LinearRoutingTable& table, std::map<std::pair<int, int> , std::vector<int>>& path_cache) {
        //distributeDemands(iter, table, path_cache);
        buildLinearRoutingFromSpanningTree(iter, table, path_cache);
    }
    void transform(TreeIteration<FlatHST>& iter, LinearRoutingTable& table, std::map<std::pair<int, int> , std::vector<int>>& path_cache) {
        //distributeDemands(iter, table, path_cache);
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
        const std::vector<int>& path,
        std::unordered_set<UndirectedEdge, UndirectedEdgeHash>& candidate_edges
    ) const {
        if (path.size() < 2) return;

        for (size_t i = 0; i + 1 < path.size(); ++i) {
            const int a = path[i];
            const int b = path[i + 1];
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

            // std::cout << "Pushing " << lambda <<" flow for commodity: " << s << std::endl;
            int cur = s;
            while (cur != root) {
                const int par = parent[cur];
                assert(par != -1);

                const int e = graph.getEdgeId(cur, par);
                // std::cout << "Edge : ( " << cur << " / " << par << " )";
                assert(e != INVALID_EDGE_ID);

                table.addFlow(e, s, lambda);
                cur = par;
            }
            // std::cout << std::endl;
        }
    }

    // ------------------------------------------------------------
    // Pointer-HST: collect all realized parent-child center paths.
    // ------------------------------------------------------------
    void collectCandidateEdges(
        const TreeIteration<std::shared_ptr<HSTNode>>& iter,
        std::unordered_set<UndirectedEdge, UndirectedEdgeHash>& candidate_edges,
        std::map<std::pair<int, int>, std::vector<int>> path_cache
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

                const auto& path = path_cache[{parentCenter, childCenter}];
                if (path.empty()) {
                    // Cache miss: compute the path and store it in the cache.
                    const auto computed_path = graph.getShortestPath(parentCenter, childCenter, distance);
                    path_cache[{parentCenter, childCenter}] = computed_path;
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
        std::map<std::pair<int, int>, std::vector<int>> path_cache
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

                const auto& path = path_cache[{parentCenter, childCenter}];
                if (path.empty()) {
                    // Cache miss: compute the path and store it in the cache.
                    const auto computed_path = graph.getShortestPath(parentCenter, childCenter, distance);
                    path_cache[{parentCenter, childCenter}] = computed_path;
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
        std::map<std::pair<int, int>, std::vector<int>>& path_cache
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

    // ------------------------------------------------------------------
    // Shared helper: apply flow along `path` (= 0 → childCenter from SPT)
    // for every (src ∈ childMembers, dst ∈ V\A) pair where src==0 or dst==0.
    // ------------------------------------------------------------------
    void applyFlow(const std::vector<int>&        path,
                   const std::unordered_set<int>&        parent_set,
                   const std::unordered_set<int>& child_member,
                   double                         lambda,
                   LinearRoutingTable&            table) {

        if (child_member.contains(0)) {
            for (int dst : parent_set) {
                if (dst == 0) continue;
                //std::cout << "Destination node: " << dst << std::endl;
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int e      = graph.getEdgeId(path[i],   path[i+1]);
                    // print flow
                    //std::cout << "Adding flow from 0 to " << dst << " along edge (" << path[i] << ", " << path[i+1] << ") with lambda " << lambda << std::endl;
                    // first check if there is already flow pushed along its reverse direction
                    int rev_edge = graph.getAntiEdge(e);
                    double rev_flow = table.getFlow(rev_edge, dst);
                    if (rev_flow > 0) {
                        //std::cout << "Existing flow in reverse direction: " << rev_flow << std::endl;
                        // subtract the flow contribution along anti edge e' by lambda
                        if (rev_flow > lambda) {
                            table.addFlow(rev_edge, dst, -lambda);
                            //std::cout << "Subtracting flow " << lambda << " from reverse edge (" << path[i+1] << ", " << path[i] << ")\n";
                        }else {
                            if (std::abs(rev_flow-lambda)<EPS) {
                                // we push the same amount of flow, so we can remove the flow contribution along e entirely
                                table.eraseFlow(rev_edge, dst);
                                //std::cout << "Removing flow entirely: " << rev_edge << " from reverse edge (" << path[i+1] << ", " << path[i] << ")\n";
                            }else {
                                // lambda is larger, so we remove flow contribution along but only push the net flow along e
                                double net_flow = lambda - rev_flow;
                                table.eraseFlow(rev_edge, dst);
                                table.addFlow(e, dst, net_flow);
                                //std::cout << "Subtracting flow " << rev_flow << " from reverse edge (" << path[i+1] << ", " << path[i] << ") and adding net flow " << net_flow << std::endl;
                            }
                        }
                    }else {
                        table.addFlow(e, dst, lambda);
                    }
                }
            }
        }

        if (parent_set.contains(0)) {
            for (int src : child_member) {
                if (src == 0) continue;
                //std::cout << "Source node: " << src << std::endl;
                for (int i = 0; i+1< path.size(); i++) {
                    int rev_e = graph.getEdgeId(path[i+1], path[i]);
                    //std::cout << "Adding flow from 0 to " << src << " along edge (" << path[i+1] << ", " << path[i] << ") with lambda " << lambda << std::endl;
                    int rev_edge = graph.getAntiEdge(rev_e);
                    double rev_flow = table.getFlow(graph.getAntiEdge(rev_e), src);
                    //table.addFlow(rev_e, src, lambda);
                    if (rev_flow > 0) {
                        //std::cout << "Existing flow in reverse direction: " << rev_flow << std::endl;
                        // subtract the flow contribution along anti edge e' by lambda
                        if (rev_flow > lambda) {
                            table.addFlow(rev_edge, src, -lambda);
                            //std::cout << "Subtracting flow " << lambda << " from reverse edge (" << path[i+1] << ", " << path[i] << ")\n";
                        }else {
                            if (std::abs(rev_flow-lambda)<EPS) {
                                // we push the same amount of flow, so we can remove the flow contribution along e entirely
                                table.eraseFlow(rev_edge, src);
                                //std::cout << "Removing flow entirely: " << rev_edge << " from reverse edge (" << path[i+1] << ", " << path[i] << ")\n";
                            }else {
                                // lambda is larger, so we remove flow contribution along but only push the net flow along e
                                double net_flow = lambda - rev_flow;
                                table.eraseFlow(rev_edge, src);
                                table.addFlow(rev_e, src, net_flow);
                                //std::cout << "Subtracting flow " << rev_flow << " from reverse edge (" << path[i+1] << ", " << path[i] << ") and adding net flow " << net_flow << std::endl;
                            }
                        }
                    }else {
                        table.addFlow(rev_e, src, lambda);
                    }
                }
            }
        }
        /*
        for (int src : childMembers) {
            for (int dst : graph.getVertices()) {
                if (A.count(dst) || src == dst) continue;
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int e      = graph.getEdgeId(path[i],   path[i+1]);
                    int anti_e = graph.getEdgeId(path[i+1], path[i]);
                    assert(e != INVALID_EDGE_ID && anti_e != INVALID_EDGE_ID);
                    if (dst == 0) {
                        table.addFlow(anti_e, src, lambda);
                    }
                    if (src == 0) {
                        table.addFlow(e,      dst, lambda);
                    }

                }
            }
        }*/
    }

    // ------------------------------------------------------------------
    // Overload 1: pointer-based HSTNode tree
    // WITH PATH CACHING to handle arbitrary parent-child center pairs
    // ------------------------------------------------------------------
    void distributeDemands(TreeIteration<std::shared_ptr<HSTNode>>& iter,
                           LinearRoutingTable& table, std::map<std::pair<int,int>, std::vector<int>>& path_cache) {
        const HSTNode*                       root       = iter.getTree().get();
        const double                         lambda     = iter.getLambda();
        const auto distance = iter.getDistance();

        int n = 0;
        for (auto v : root->getMembers()) n = std::max(n, v);
        n++;

        // --- Cache for computed shortest paths (parentCenter → childCenter) ---
        // Use map to avoid hash collisions; key = (parentCenter, childCenter)
        //std::map<std::pair<int,int>, std::vector<int>> path_cache;

        std::vector<bool> child_membership(n, false);

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

                // --- Check cache first, compute and cache if not found ---
                auto cache_key = std::make_pair(parentCenter, childCenter);
                auto& path = path_cache[cache_key];
                if (path.empty()) {
                    path = graph.getShortestPath(parentCenter, childCenter, distance);
                }

                if (path.size() < 2) { q.push(child); continue; }

                // --- Build membership vector for O(1) lookups ---
                std::fill(child_membership.begin(), child_membership.end(), false);
                for (int member : child->getMembers()) {
                    child_membership[member] = true;
                }

                applyFlowOptimized(path, child_membership, lambda, table, n);
                q.push(child);
            }
        }
    }

    // ------------------------------------------------------------------
    // Overload 2: flat FlatHST tree
    // WITH PATH CACHING to handle arbitrary parent-child center pairs
    // ------------------------------------------------------------------
    void distributeDemands(TreeIteration<FlatHST>& iter,
                           LinearRoutingTable& table,
                           std::map<std::pair<int,int>, std::vector<int>>& path_cache) {
        const FlatHST&                           hst        = iter.getTree();
        const double                         lambda     = iter.getLambda();
        const auto distance = iter.getDistance();

        int n = hst.nodes.size();

        // --- Cache for computed shortest paths (parentCenter → childCenter) ---
        // Use map to avoid hash collisions; key = (parentCenter, childCenter)


        std::vector<bool> child_membership(n, false);

        //print(hst);
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

                // --- Check cache first, compute and cache if not found ---
                auto cache_key = std::make_pair(parentCenter, childCenter);
                auto& path = path_cache[cache_key];
                if (path.empty()) {
                    path = graph.getShortestPath(parentCenter, childCenter, distance);
                }

                if (path.size() < 2) { q.push(child_idx); continue; }

                // --- Build membership vector for O(1) lookups ---
                std::fill(child_membership.begin(), child_membership.end(), false);
                for (int member : child_members) {
                    child_membership[member] = true;
                }

                applyFlowOptimized(path, child_membership, lambda, table, n);
                q.push(child_idx);
            }
        }
    }

    // ------------------------------------------------------------------
    // Optimized flow application using vector<bool> for O(1) membership tests
    // and precomputed edge IDs to avoid redundant lookups
    // ------------------------------------------------------------------
    void applyFlowOptimized(const std::vector<int>&  path,
                            const std::vector<bool>& child_membership,
                            double                   lambda,
                            LinearRoutingTable&      table,
                            int                      n) {

        // --- Precompute all edge IDs on the path to avoid repeated lookups ---
        std::vector<int> edge_ids;
        std::vector<int> rev_edge_ids;
        edge_ids.reserve(path.size() - 1);
        rev_edge_ids.reserve(path.size() - 1);

        for (size_t i = 0; i + 1 < path.size(); ++i) {
            int e = graph.getEdgeId(path[i], path[i+1]);
            int rev_e = graph.getEdgeId(path[i+1], path[i]);
            edge_ids.push_back(e);
            rev_edge_ids.push_back(rev_e);
        }

        bool child_has_0 = child_membership[0];

        // --- Case 1: Child contains node 0 → push flow towards other nodes ---
        if (child_has_0) {
            for (int dst = 0; dst < n; ++dst) {
                if (dst == 0 || child_membership[dst]) continue;

                for (size_t i = 0; i < edge_ids.size(); ++i) {
                    int e = edge_ids[i];
                    int rev_edge = graph.getAntiEdge(e);
                    double rev_flow = table.getFlow(rev_edge, dst);

                    if (rev_flow > EPS) {
                        if (rev_flow > lambda + EPS) {
                            table.addFlow(rev_edge, dst, -lambda);
                        } else if (std::abs(rev_flow - lambda) < EPS) {
                            table.eraseFlow(rev_edge, dst);
                        } else {
                            double net_flow = lambda - rev_flow;
                            table.eraseFlow(rev_edge, dst);
                            table.addFlow(e, dst, net_flow);
                        }
                    } else {
                        table.addFlow(e, dst, lambda);
                    }
                }
            }
        }

        // --- Case 2: Child does NOT contain 0 → push flow from child nodes to root ---
        if (!child_has_0) {
            for (int src = 0; src < n; ++src) {
                if (src == 0 || !child_membership[src]) continue;

                for (size_t i = 0; i < rev_edge_ids.size(); ++i) {
                    int rev_e = rev_edge_ids[i];
                    int rev_edge = graph.getAntiEdge(rev_e);
                    double rev_flow = table.getFlow(rev_edge, src);

                    if (rev_flow > EPS) {
                        if (rev_flow > lambda + EPS) {
                            table.addFlow(rev_edge, src, -lambda);
                        } else if (std::abs(rev_flow - lambda) < EPS) {
                            table.eraseFlow(rev_edge, src);
                        } else {
                            double net_flow = lambda - rev_flow;
                            table.eraseFlow(rev_edge, src);
                            table.addFlow(rev_e, src, net_flow);
                        }
                    } else {
                        table.addFlow(rev_e, src, lambda);
                    }
                }
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
