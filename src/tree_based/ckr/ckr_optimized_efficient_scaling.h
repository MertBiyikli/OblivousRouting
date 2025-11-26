//
// Created by Mert Biyikli on 03.11.25.
//

#ifndef OBLIVIOUSROUTING_CKR_OPTIMIZED_EFFICIENT_SCALING_H
#define OBLIVIOUSROUTING_CKR_OPTIMIZED_EFFICIENT_SCALING_H

#include "../../solver/solver.h"
#include <random>
#include "ckr_tree_decomposer.h"
#include "utils/ultrametric_tree.h"
#include "raecke_ckr_transform.h"

#include <chrono>

/*
 * Optimized Raecke CKR implementation with efficient scaling:
 * - use tighter scaling factors
 * - use the precomputed MST as preprocessing
 * - embrace LCA DS for better performance and Quotient Graph Generating
 */
namespace MendelScaling {


    class CKRTransform;
    class RaeckeCKR;



    /*
     * Implement the RaeckeCKR_EfficientScaling
     */
    class RaeckeCKR {
        // similar members as RaeckeCKR but with efficient scaling modifications
    public:


        void run() {
            m_lambdaSum = 0.0;
            int id = 0;
            while (m_lambdaSum < 1.0) {
                auto start = std::chrono::high_resolution_clock::now();
                m_lambdaSum += iterate(id);

                oracle_running_times.push_back(
                    std::chrono::duration<double, std::milli>(
                        std::chrono::high_resolution_clock::now() - start
                    ).count()
                );
                id++;
            }
        }


        double iterate(int id) {
            std::shared_ptr<TreeNode> t = getTree(m_graph);
            computeRLoads(t, m_graph, id);
            double l = getMaxRload(id);
            double lambda = std::min(1.0/l, 1.0 - m_lambdaSum);


            // for debuggin purposes print the tree , graph and the current edge loads, as well as the current lambda
            if (debug) {
                std::cout << "Iteration " << id << ":\n";
                std::cout << "Current lambda: " << lambda << "\n";
                std::cout << "Current Tree:\n";
                print_tree(t);
                std::cout << "Current Graph Distances:\n";
                m_graph.print();
                std::cout << "Current Edge R-Loads:\n";
                for (const auto& [edge, rLoad] : edge2Load[id]) {
                    std::cout << "Edge (" << edge.first << ", " << edge.second << ") : R-Load = " << rLoad << "\n";
                }
            }


            m_lambdas.push_back(lambda);
            m_graphs.push_back(m_graph);
            m_trees.push_back(t);
            // update weights
            m_lambdaSum += lambda;
            computeNewDistances(m_graph);

            edge2Load.clear(); // clear for next iteration
            return lambda;
        }

        const int getIterationCount() const {
            return static_cast<int>(oracle_running_times.size());
        }


        std::vector<CKRLevel> m_levels;          // from finest (0) upward


        std::vector<int> build_ckr_level(const Graph& g, double Delta, CKRLevel& L);

        /*
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
        */

    bool debug = false;

        Graph m_graph;
        double m_lambdaSum = 0.0;
        double diameter = 0.0;


        // TODO: keep this for now, but can lead to memory performance issue
        std::vector<Graph> m_graphs;
        std::vector<double> m_lambdas;
        std::vector<std::shared_ptr<TreeNode>> m_trees;
        UltrametricTree ultrametric;

        bool use_mendel_scaling = true;      // toggle at runtime or via CLI
        uint32_t ckr_seed = 0xC0FFEE;       // for reproducible partitions


        std::vector<std::unordered_map<std::pair<int, int>, double>> edge2Load;
        std::vector<double> oracle_running_times;
        void init(const Graph& g);
        void preprocess();
        std::shared_ptr<TreeNode> getTree(Graph& g);
        void computeRLoads(std::shared_ptr<TreeNode> t, Graph& g, int tree_index);
        double getMaxRload(int tree) const;
        // void addLoadToEdge(int u, int v, double load);
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
        void addTree(std::shared_ptr<TreeNode> root, double lambda, Graph& graph);
    };

    class RaeckeCKROptimized : public ObliviousRoutingSolver
    {


    public:
        RaeckeCKR ckr_algo;
        RaeckeCKRTransform transform;
        RaeckeCKROptimized() = default;
        ~RaeckeCKROptimized() = default;

        void solve(const Graph &g) override {
            ckr_algo.debug = debug;
            ckr_algo.setGraph(g);

            ckr_algo.init(g);
            ckr_algo.run();  // pass transform directly
            iteration_count = ckr_algo.getIterationCount();
            oracle_running_times = ckr_algo.oracle_running_times;
            scaleDownFlow(); // normalize after building
        }

        void storeFlow() override {
            // directly add flow contribution
            for (int i = 0; i < ckr_algo.m_graphs.size(); ++i) {
                transform.addTree(ckr_algo.m_trees[i],ckr_algo.m_lambdas[i], ckr_algo.m_graphs[i]);
            }

            // store the flow
            // given the demand map
            auto const& routingRaecke = transform.getRoutingTable();
            for (const auto& [edge, demandMap] : routingRaecke) {
                for (const auto& [d, fraction] : demandMap) {
                    if (d.first > d.second) continue; // skip trivial cases
                    f_e_st[edge][d]=fraction;
                }
            }

            scaleDownFlow();
/*
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
            }*/
        } // handled during run()

        void scaleDownFlow() {
            // scale the flow to meet unit flow
            std::unordered_map<std::pair<int, int>, double > outgoingflow_per_commodity;

            for ( const auto& [edge, flowMap]:f_e_st) {
                for (const auto& [com, flow_value]  : flowMap) {
                    if ( flow_value < 1e-15 ) continue; // ignore zero flows
                    if (!outgoingflow_per_commodity.contains(com) )
                        outgoingflow_per_commodity[com] = 0;

                    if (edge.first == com.first
                        || edge.second == com.first) {
                        outgoingflow_per_commodity[com] += std::abs(flow_value);
                        }
                }
            }

            // scale the flow values to meet one unit of flow per commodity
            for ( auto& [edge, flowMap]:f_e_st) {
                for (auto& [com, flow_value] : flowMap) {
                    flow_value /= outgoingflow_per_commodity[com];
                }
            }
        }
    };
}


#endif //OBLIVIOUSROUTING_CKR_OPTIMIZED_EFFICIENT_SCALING_H