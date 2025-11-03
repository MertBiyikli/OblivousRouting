//
// Created by Mert Biyikli on 27.10.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_CKR_OPTIMIZED_H
#define OBLIVIOUSROUTING_RAECKE_CKR_OPTIMIZED_H


#include <random>

#include "../../solver/solver.h"
#include "ckr_tree_decomposer.h"

class RaeckeCKR {
    Graph m_graph;
    double m_lambdaSum = 0.0;
public:
    // store the current load of the tree for the edges+
    /*
    std::vector<std::vector<int>> m_adj_edge_loads;
    std::vector<std::vector<double>> m_edgeLoads;*/

    std::unordered_map<std::pair<int, int>, double> edge2Load;

    std::vector<double> oracle_running_times;
    void init(Graph& g);
    TreeNode* getTree(Graph& g);
    void computeRLoads(int idx, TreeNode* t, Graph& g);
    double getMaxRload();
    void addLoadToEdge(int u, int v, double load);
    void computeNewDistances(Graph& g);
    void setGraph(const Graph& g) {
        m_graph = g;
    }

    const Graph& getGraph() const {
        return m_graph;
    }

    template<class T>
    void run(T& transform) {
        m_lambdaSum = 0.0;
        int id = 0;
        while (m_lambdaSum < 1.0) {
            auto start = std::chrono::high_resolution_clock::now();
            m_lambdaSum += iterate(id, transform);

            oracle_running_times.push_back(
                std::chrono::duration<double, std::milli>(
                    std::chrono::high_resolution_clock::now() - start
                ).count()
            );
            id++;
        }
    }

    template<class T>
    double iterate(int id, T& transform) {
        TreeNode* t = getTree(m_graph);
        computeRLoads(id, t, m_graph);
        double l = getMaxRload();
        double lambda = std::min(1.0/l, 1.0 - m_lambdaSum);

        // directly add flow contribution
        transform.addTree(t, lambda, m_graph);

        // update weights
       // m_lambdaSum += lambda;
        computeNewDistances(m_graph);

        edge2Load.clear(); // clear for next iteration
        return lambda;
    }

    const int getIterationCount() const {
        return static_cast<int>(oracle_running_times.size());
    }

    struct CKRLevel {
        double R = 0.0;                      // the random radius used at this level
        std::vector<int> owner;              // owner[v] = center that captured v at this level (or -1)
        std::vector<int> pred;               // pred[v] = predecessor of v towards its owner center at this level (-1 for center)
        std::vector<int> centers;            // the centers chosen at this level (subset of V)
    };

    std::vector<CKRLevel> m_levels;          // from finest (0) upward

    // Build one CKR level at scale Delta. Returns next-level centers.
static std::vector<int> build_ckr_level(
    const Graph& G,
    double Delta,
    CKRLevel& L,
    std::mt19937& rng
) {
    const int n = G.getNumNodes();
    L.owner.assign(n, -1);
    L.pred.assign(n, -1);

    // R ~ Uniform[Delta/4, Delta/2]
    std::uniform_real_distribution<double> U(0.25 * Delta, 0.5 * Delta);
    L.R = U(rng);

    // Random permutation of all vertices as candidate centers
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);

    // Truncated Dijkstra from centers in perm order; skip already owned verts
    using P = std::pair<double,int>;
    std::priority_queue<P, std::vector<P>, std::greater<P>> pq;

    std::vector<double> dist(n, std::numeric_limits<double>::infinity());
    std::vector<int>    dist_epoch(n, 0);
    int epoch = 1;

    auto reset_heap_and_epoch = [&]() {
        while (!pq.empty()) pq.pop();
        ++epoch;  // lazy-reset dist[] via epoch
    };
    auto has_dist = [&](int v) { return dist_epoch[v] == epoch; };
    auto set_dist = [&](int v, double d) { dist_epoch[v] = epoch; dist[v] = d; };

    L.centers.clear();
    L.centers.reserve(n);

    for (int c : perm) {
        if (L.owner[c] != -1) continue;          // already captured
        // Start a new center
        L.centers.push_back(c);
        reset_heap_and_epoch();
        set_dist(c, 0.0);
        L.pred[c] = -1;
        pq.emplace(0.0, c);

        while (!pq.empty()) {
            auto [d,u] = pq.top(); pq.pop();
            if (d >= L.R) break;
            if (!has_dist(u) || d != dist[u]) continue;
            if (L.owner[u] != -1) continue;      // skip: already assigned by earlier center

            // Capture u for center c
            L.owner[u] = c;

            for (auto& v : G.neighbors(u)) {
                if (L.owner[v] != -1) continue;  // already captured
                double w = G.getEdgeDistance(u, v);
                double nd = d + w;
                if (nd < L.R && (!has_dist(v) || nd < dist[v])) {
                    set_dist(v, nd);
                    L.pred[v] = u;               // predecessor towards center c
                    pq.emplace(nd, v);
                }
            }
        }
    }
    return L.centers; // these are the "surviving" centers for the next (coarser) scale
}
};


struct FlowStorage {

    int m = 0, n = 0;
    void init(const Graph& g) {
        m = g.getNumEdges();
        n = g.getNumNodes();
        adj_f_e_u_id.resize(m);
        adj_f_e_u.resize(m);

        // init the edges list
        edges.reserve(m);
        for (int v = 0; v<g.getNumNodes(); v++) {
            for (const auto& u : g.neighbors(v)) {
                if (v < u) {
                    edges.emplace_back(v,u);
                }
            }
        }

        // sort the edges
        std::sort(edges.begin(), edges.end(),
                  [](const std::pair<int,int> &a, const std::pair<int,int> &b) {
                      if (a.first != b.first) return a.first < b.first;
                      return a.second < b.second;
                  });
    }

    void addFlow(int u, int v, int s, double flow) {
        // call the optimized addFlow
        int a = std::min(u,v);
        int b = std::max(u,v);

        // find edge id
        int lo = 0, hi = m;
        while (lo < hi) {
            int mid = (lo + hi) >> 1;
            const auto& [eu, ev] = edges[mid];
            if (eu < a || (eu == a && ev < b))
                lo = mid + 1;
            else
                hi = mid;
        }

        // call optimized addFlow
        // if (lo < m && edges[lo] == std::make_pair(a,b)) {
            addFlow(lo, s, flow);
        // } else {
            // edge not found (should not happen)
           //  throw std::runtime_error("Edge not found in FlowStorage::addFlow");
        // }
    }

// TODO: make the adjacency list representation global for all solvers, except maybe the Applegate and Cohen one...
    std::vector<std::pair<int,int>> edges; // edge list (u,v) with u < v
    std::vector<std::vector<int>> adj_f_e_u_id; // per-adjacency list version of f_e_u (stores edge ids)
    std::vector<std::vector<double>> adj_f_e_u; // per-adjacency list version of f_e_u
    // -----------------------------------------------------------------------------
    // ⚡ Optimized addFlow (manual binary search, cache-friendly)
    // -----------------------------------------------------------------------------
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((always_inline))
    #elif defined(_MSC_VER)
    __forceinline
    #endif
    inline void addFlow(int edge_id, int u, double flow) noexcept {
        auto &ids  = adj_f_e_u_id[edge_id];
        auto &vals = adj_f_e_u[edge_id];
        const size_t len = ids.size();

        // Reserve to avoid reallocations
        if (len == 0) [[unlikely]] {
            ids.reserve(8);
            vals.reserve(8);
            ids.push_back(u);
            vals.push_back(flow);
            return;
        }

        // ✅ Linear scan for small adjacency lists (<= 16 entries)
        if (len <= 16) [[likely]] {
            for (size_t i = 0; i < len; ++i) {
                if (ids[i] == u) {
                    vals[i] += flow;
                    return;
                }
            }
            ids.push_back(u);
            vals.push_back(flow);
            return;
        }

        // ✅ Hand-rolled binary search (faster than std::lower_bound)
        size_t lo = 0, hi = len;
        while (lo < hi) {
            const size_t mid = (lo + hi) >> 1;
            const int mid_val = ids[mid];
            if (mid_val < u)
                lo = mid + 1;
            else
                hi = mid;
        }

        if (lo < len && ids[lo] == u) {
            vals[lo] += flow;
        } else {
            // Usually append to the end — amortized O(1)
            if (lo == len) {
                ids.push_back(u);
                vals.push_back(flow);
            } else {
                ids.insert(ids.begin() + static_cast<long>(lo), u);
                vals.insert(vals.begin() + static_cast<long>(lo), flow);
            }
        }
    }

    void scaleDownFlow() {
        // do the recaling for the adjacency matrix as well
        std::vector<double> outgoingflow_f_s_x_fixed(n, 0);
        for (int e = 0; e < m; ++e) {
            const int u = edges[e].first;
            const int v = edges[e].second;
            const auto &ids  = adj_f_e_u_id[e];
            const auto &vals = adj_f_e_u[e];
            const int L = static_cast<int>(ids.size());

            for (int k = 0; k < L; ++k) {
                int s = ids[k];
                double f = vals[k];

                if ( f >= 0 ) {
                    // direction u -> v
                    if (s == u) {
                        outgoingflow_f_s_x_fixed[s] += f;
                    }
                } else {
                    // direction v -> u
                    if (s == v) {
                        outgoingflow_f_s_x_fixed[s] += std::abs(f);
                    }
                }
            }
        }


        // ---------- 3) Rescale all (e, s) by the commodity's net outflow ----------
        // After this, each commodity s should have unit net outflow at s.
        constexpr double EPS = 1e-12;
        for (int e = 0; e < m; ++e) {
            auto &ids  = adj_f_e_u_id[e];
            auto &vals = adj_f_e_u[e];
            const int L = static_cast<int>(ids.size());

            for (int k = 0; k < L; ++k) {
                const int s = ids[k];

                const double denom = outgoingflow_f_s_x_fixed[s];
                if (std::abs(denom) <= EPS) {
                    if ( true) {
                        std::cerr << "[scaleFlowDown] Warning: net_outflow ~ 0 for source "
                                  << s << " (|denom|=" << std::abs(denom) << "), skipping rescale.\n";
                    }
                    continue;
                }

                const double scaled = vals[k] / denom;
                vals[k] = scaled;

                // Keep the old map in sync for downstream code (e.g., getFlowForCommodity).
                // The map key is (edge_id, source).
                //f_e_u[{e, s}] = scaled;
            }
        }
    }
};

class CKRTransform {
public:

    FlowStorage flow_storage;
    // edge_to_id / id_to_edge are no longer needed unless you want edge indexing later
    std::unordered_map<std::pair<int, int>, std::unordered_map<std::pair<int,int>, double>> flow_map;

    void init(const Graph& g) {
        // fix a node for the storing the lienar oblivious routing
        this->flow_storage.init(g);
    }

    inline void addFlow(int u, int v, int s, int t, double flow) {
        flow_storage.addFlow(u, v, s, flow);
        // flow_map[{u, v}][{s, t}] += std::abs(flow);
        // int a = std::min(u,v), b = std::max(u,v);
    }

    void addTree(const TreeNode* root, double lambda, Graph& graph) {
        if (!root) return;
        std::queue<const TreeNode*> q; q.push(root);

        while (!q.empty()) {
            const auto* node = q.front(); q.pop();
            for (const auto* child : node->children) {
                q.push(child);
                const auto& S = child->members;
                if (S.empty()) continue;

                // --- Precompute membership mask for S ---
                std::vector<char> inS(graph.getNumNodes(), 0);
                for (int v : S) inS[v] = 1;

                int repParent = node->members.empty() ? S[0] : node->members[0];
                int repChild  = S[0];

                // TODO: improve representative selection
                //  ideally use something like an LCA data structure to avoid shortest path computations
                auto path = graph.getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                const double per_source = lambda;// / (double)S.size();

                // Each vertex in S sends λ flow to the fixed node
                for (int s : S) {
                    for (int t = 0; t < graph.getNumNodes(); ++t) {
                        if (inS[t] || s == t) continue;
                        for (size_t i = 0; i + 1 < path.size(); ++i) {
                            int a = path[i];
                            int b = path[i + 1];
                            double sign = (a < b) ? +1.0 : -1.0;


                            addFlow(a, b, s, t, sign * per_source);
                        }
                    }
                }
            }
        }
    }
};

class RaeckeCKROptimized : public ObliviousRoutingSolver {
    RaeckeCKR ckr_algo;
    CKRTransform transform;

public:

    RaeckeCKROptimized() = default;
    ~RaeckeCKROptimized() = default;

    void solve(const Graph &g) override {
        ckr_algo.setGraph(g);
        transform.init(g);
        ckr_algo.run(transform);  // pass transform directly
        iteration_count = ckr_algo.getIterationCount();
        oracle_running_times = ckr_algo.oracle_running_times;
        scaleDownFlow(); // normalize after building
    }

    void storeFlow() override {
        f_e_st.clear();
        // store the flow from the adjacency list flow
        for (int e = 0; e<transform.flow_storage.adj_f_e_u_id.size(); e++) {
            int u = transform.flow_storage.edges[e].first;
            int v = transform.flow_storage.edges[e].second;
            for (int s = 0; s<ckr_algo.getGraph().getNumNodes(); ++s) {
                for (int t = 0; t<ckr_algo.getGraph().getNumNodes(); t++) {
                    if ( s == t ) continue;
                    double flow(0);
                    auto it = std::find(transform.flow_storage.adj_f_e_u_id[e].begin(),
                        transform.flow_storage.adj_f_e_u_id[e].end(), s);
                    if (it != transform.flow_storage.adj_f_e_u_id[e].end()) {
                        int idx = std::distance(transform.flow_storage.adj_f_e_u_id[e].begin(), it);
                        flow = transform.flow_storage.adj_f_e_u[e][idx];
                    }
                    else {
                        auto it2 = std::find(transform.flow_storage.adj_f_e_u_id[e].begin(),
                            transform.flow_storage.adj_f_e_u_id[e].end(), t);
                        if (it2 != transform.flow_storage.adj_f_e_u_id[e].end()) {
                            int idx = std::distance(transform.flow_storage.adj_f_e_u_id[e].begin(), it2);
                            flow = -transform.flow_storage.adj_f_e_u[e][idx];
                        }
                    }
                    if (std::abs(flow) > 1e-12) {
                        // since this linear
                        if (flow <0) {
                            f_e_st[{v, u}][{s, t}] = -flow; // Store the reverse flow as well
                        }else {
                            f_e_st[{u, v}][{s, t}] = flow;
                        }
                    }
                }
            }
        }
    } // handled during run()

    void scaleDownFlow() {
        transform.flow_storage.scaleDownFlow();

/*
        std::unordered_map<std::pair<int,int>, double> total_outflow;
        auto routing = transform.flow_map;
        for (const auto& [edge, dmap] : routing) {
            for (const auto& [dem, flow] : dmap) {
                if (edge.first == dem.first ||
                    edge.second == dem.first)
                    total_outflow[dem] += std::abs(flow);
            }
        }

        // print out the outgoing flows
        for (const auto& [dem, total_flow] : total_outflow) {
            if (total_flow < 1e-12) {
                std::cout << "Warning: Total outflow for commodity (" << dem.first << ", " << dem.second << ") is zero or very small." << std::endl;
            }
        }
        for (auto& [edge, dmap] : routing) {
            for (auto& [dem, flow] : dmap)
                routing[edge][dem] = flow / total_outflow[dem];
        }
        f_e_st = routing; // final flow table
*/
    }
};

#endif //OBLIVIOUSROUTING_RAECKE_CKR_OPTIMIZED_H