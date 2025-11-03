//
// Created by Mert Biyikli on 03.11.25.
//

#ifndef OBLIVIOUSROUTING_CKR_OPTIMIZED_EFFICIENT_SCALING_H
#define OBLIVIOUSROUTING_CKR_OPTIMIZED_EFFICIENT_SCALING_H

#include "../../solver/solver.h"
#include <random>
#include "ckr_tree_decomposer.h"

/*
 * Optimized Raecke CKR implementation with efficient scaling:
 * - use tighter scaling factors
 * - use the precomputed MST as preprocessing
 * - embrace LCA DS for better performance and Quotient Graph Generating
 */
namespace MendelScaling {


    class CKRTransform;
    class RaeckeCKR;
    struct CKRLevel {
        double R = 0.0;                      // the random radius used at this level
        std::vector<int> owner;              // owner[v] = center that captured v at this level (or -1)
        std::vector<int> pred;               // pred[v] = predecessor of v towards its owner center at this level (-1 for center)
        std::vector<int> centers;            // the centers chosen at this level (subset of V)
    };

struct UltrametricTree {
    // Each original vertex v is a leaf at index v (0..n-1).
    // Internal nodes are appended at indices n .. n+merges-1.
    struct Node {
        int parent = -1;
        std::vector<int> child;
        double Gamma = 0.0;   // label Γ(u) (nondecreasing toward root)
    };

    int n = 0;                // #original vertices
    int N = 0;                // total nodes = leaves + internal
    int root = -1;
    std::vector<Node> T;

    // Binary lifting tables for level-ancestor-like jumps
    std::vector<std::vector<int>> up;     // up[k][u] = 2^k-th ancestor of u (or -1)
    std::vector<double> gamma;            // Γ(u), cached for easy access
    std::vector<int> depth;               // depth in the ultrametric tree

    // Build ultrametric from an MST using Kruskal-like dendrogram construction:
    // Sort MST edges by weight ascending; each union creates a new internal node
    // with Γ = edge weight; connect components as its children.
    template<class EdgeList>
    void buildFromMST(int n_, const EdgeList& mst_edges, const std::vector<double>& mst_w) {
        n = n_;
        const int m = (int)mst_edges.size();
        // DSU
        std::vector<int> parent(2*n + m, -1), dsu(2*n + m);
        std::iota(dsu.begin(), dsu.end(), 0);
        auto find = [&](int x) {
            while (dsu[x] != x) { dsu[x] = dsu[dsu[x]]; x = dsu[x]; }
            return x;
        };
        auto unite = [&](int a, int b) {
            a = find(a); b = find(b);
            if (a == b) return -1;
            dsu[a] = a; dsu[b] = b; // make sure reps are clean
            return 0;
        };

        // Initially each leaf is its own component represented by itself (0..n-1)
        // We will create internal nodes starting at id = cur = n
        T.clear(); T.resize(n); // n leaf nodes
        for (int v = 0; v < n; ++v) T[v].Gamma = 0.0;

        // We need a separate DSU parent array that points to current representative node id.
        // rep[v] maps DSU representative to the current *tree node id* of that component.
        std::vector<int> rep(2*n + m);
        for (int v = 0; v < n; ++v) rep[v] = v, dsu[v] = v;

        // Sort MST edges by weight ascending
        std::vector<int> ord(m);
        std::iota(ord.begin(), ord.end(), 0);
        std::sort(ord.begin(), ord.end(), [&](int i, int j){ return mst_w[i] < mst_w[j]; });

        int cur = n; // next internal node id
        for (int idx : ord) {
            int u = mst_edges[idx].first;
            int v = mst_edges[idx].second;
            double w = mst_w[idx];

            int ru = find(rep[u]);
            int rv = find(rep[v]);
            if (ru == rv) continue;

            // Create new internal node
            if ((int)T.size() <= cur) T.resize(cur+1);
            T[cur].Gamma = w;
            int nodeU = rep[ru];
            int nodeV = rep[rv];
            T[cur].child.push_back(nodeU);
            T[cur].child.push_back(nodeV);
            T[nodeU].parent = cur;
            T[nodeV].parent = cur;

            // Unite DSU reps and map to new node
            dsu[ru] = cur;
            dsu[rv] = cur;
            rep[cur] = cur;
            rep[ru] = cur;
            rep[rv] = cur;

            ++cur;
        }
        N = cur;
        // Find root (the last internal or the single leaf)
        root = (N > 0 ? N-1 : (n>0? n-1 : -1));
        if (root == -1 && n > 0) root = n-1;

        // Preprocess lifting
        preprocessLifting();
    }

    void preprocessLifting() {
        gamma.assign(N, 0.0);
        for (int i = 0; i < N; ++i) gamma[i] = T[i].Gamma;

        depth.assign(N, 0);
        std::queue<int> q;
        if (root != -1) q.push(root), depth[root] = 0;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : T[u].child) {
                depth[v] = depth[u] + 1;
                q.push(v);
            }
        }
        int LOG = 1;
        while ((1<<LOG) <= std::max(1, N)) ++LOG;
        up.assign(LOG, std::vector<int>(N, -1));
        for (int u = 0; u < N; ++u) {
            up[0][u] = T[u].parent;
        }
        for (int k = 1; k < LOG; ++k) {
            for (int u = 0; u < N; ++u) {
                int mid = up[k-1][u];
                up[k][u] = (mid == -1 ? -1 : up[k-1][mid]);
            }
        }
    }

    // Return the *highest* ancestor of leaf v whose Γ <= threshold.
    // If leaf’s Γ(leaf)=0 is already > threshold (never happens), it returns the leaf itself.
    int sigmaDelta(int v_leaf, double Delta) const {
        if (N == 0 || v_leaf < 0 || v_leaf >= n) return v_leaf;
        const double thr = Delta / (2.0 * (double) n);

        int u = v_leaf;
        if (gamma[u] > thr) return u; // safety
        // climb as long as parent exists and Γ(parent) <= thr
        for (int k = (int)up.size()-1; k >= 0; --k) {
            int p = up[k][u];
            if (p != -1 && gamma[p] <= thr) {
                u = p;
            }
        }
        return u;
    }
};

    struct QuotientLevel {
    Graph Gq;                                 // quotient graph at Δ
    std::vector<int> sigma_compact_of_v;      // size n: original vertex v -> compact qid in [0..k-1]
    std::vector<std::vector<int>> members_of_q; // size k: list of original vertices in each quotient node
};

// Build the Δ-level quotient graph and all mappings needed to map back to original vertices.
// Complexity: O(m log n) amortized over all levels (by the paper’s edge-coverage argument).
inline QuotientLevel build_quotient_graph_with_map(const Graph& G, const UltrametricTree& ultra, double Delta) {
    const int n = G.getNumNodes();

    // 1) Map each original vertex to ancestor σΔ(v)
    std::vector<int> sigma_node(n);
    for (int v = 0; v < n; ++v) {
        sigma_node[v] = ultra.sigmaDelta(v, Delta); // ultra node id
    }

    // 2) Compact those ancestor IDs to 0..k-1
    std::unordered_map<int,int> idmap;
    idmap.reserve(n*2);
    int next_id = 0;
    std::vector<int> qid_of_sigma(ultra.N, -1);
    for (int v = 0; v < n; ++v) {
        int s = sigma_node[v];
        int &ref = qid_of_sigma[s];
        if (ref == -1) ref = next_id++;
    }

    // 3) Original v -> compact quotient id
    std::vector<int> sigma_compact_of_v(n);
    sigma_compact_of_v.reserve(n);
    for (int v = 0; v < n; ++v) {
        sigma_compact_of_v[v] = qid_of_sigma[sigma_node[v]];
    }
    const int k = next_id;

    // 4) Members of each quotient node
    std::vector<std::vector<int>> members_of_q(k);
    for (int v = 0; v < n; ++v) {
        members_of_q[sigma_compact_of_v[v]].push_back(v);
    }

    // 5) Build quotient edges with min inter-edge weight
    Graph Gq(k);
    std::unordered_map<long long,double> min_w;
    min_w.reserve((size_t)G.getNumEdges());
    auto key = [](int a,int b){ return ((long long)a<<32) | (unsigned)b; };

    for (int u = 0; u < n; ++u) {
        int cu = sigma_compact_of_v[u];
        for (int v : G.neighbors(u)) {
            int cv = sigma_compact_of_v[v];
            if (cu == cv) continue;
            int a = (cu < cv) ? cu : cv;
            int b = (cu < cv) ? cv : cu;
            double w = G.getEdgeDistance(u, v);
            long long K = key(a,b);
            auto it = min_w.find(K);
            if (it == min_w.end() || w < it->second)
                min_w[K] = w;
        }
    }
    for (auto& [K, w] : min_w) {
        int a = (int)(K >> 32);
        int b = (int)(K & 0xffffffff);
        Gq.addEdge(a,b,w);
    }

    return {.Gq = std::move(Gq),
            .sigma_compact_of_v = std::move(sigma_compact_of_v),
            .members_of_q = std::move(members_of_q)};
}




    /*
     * Implement the RaeckeCKR_EfficientScaling
     */
    class RaeckeCKR {
        // similar members as RaeckeCKR but with efficient scaling modifications
    public:

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
            computeRLoads(t, m_graph);
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

    bool debug = false;

        Graph m_graph;
        double m_lambdaSum = 0.0;
        double diameter = 0.0;
        UltrametricTree ultrametric;

        bool use_mendel_scaling = true;      // toggle at runtime or via CLI
        uint32_t ckr_seed = 0xC0FFEE;       // for reproducible partitions


        std::unordered_map<std::pair<int, int>, double> edge2Load;
        std::vector<double> oracle_running_times;
        void init(const Graph& g);
        void preprocess();
        TreeNode* getTree(Graph& g);
        void computeRLoads(TreeNode* t, Graph& g);
        double getMaxRload();
        void addLoadToEdge(int u, int v, double load);
        void computeNewDistances(Graph& g);
        void setGraph(const Graph& g);
        const Graph& getGraph() const;
    };

    struct FlowStorage {

        int m = 0, n = 0;
        void init(const Graph& g);
        void addFlow(int u, int v, int s, double flow);

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
        inline void addFlow(int edge_id, int u, double flow) noexcept;

        void scaleDownFlow();
    };
    class CKRTransform {
        public:
        FlowStorage flow_storage;
        std::unordered_map<std::pair<int, int>, std::unordered_map<std::pair<int,int>, double>> flow_map;

        void init(const Graph& g);
        inline void addFlow(int u, int v, int s, int t, double flow);
        void addTree(const TreeNode* root, double lambda, Graph& graph);
    };

    class RaeckeCKROptimized : public ObliviousRoutingSolver
    {
        RaeckeCKR ckr_algo;
        CKRTransform transform;

    public:

        RaeckeCKROptimized() = default;
        ~RaeckeCKROptimized() = default;

        void solve(const Graph &g) override {
            ckr_algo.setGraph(g);
            transform.init(g);
            ckr_algo.init(g);
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
        }
    };
}


#endif //OBLIVIOUSROUTING_CKR_OPTIMIZED_EFFICIENT_SCALING_H