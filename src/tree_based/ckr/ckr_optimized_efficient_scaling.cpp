//
// Created by Mert Biyikli on 03.11.25.
//

#include "ckr_optimized_efficient_scaling.h"
#include "../random_mst/mst.h"
#include <cassert>
#include <chrono>


namespace MendelScaling {
    void RaeckeCKR::init(const Graph &g) {
        this->m_graph = g;
        for (int i = 0; i < g.getNumNodes(); ++i) {
            for (auto& u : g.neighbors(i)) {
                if ( i < u) {
                    double cap = g.getEdgeCapacity(i, u);
                    m_graph.updateEdgeDistance(i, u, cap);
                }
            }
        }
        this->m_graph.createDistanceMatrix();
        this->diameter = m_graph.GetDiameter();
    }


    void RaeckeCKR::preprocess() {

        // -- compute MST using Kruskals Algo
        RandomMST mst_oracle;
        mst_oracle.setGraph(m_graph);
        std::vector<std::pair<int, int> > mst = mst_oracle.build_mst(m_graph);


        // ----  Collect MST edge weights ----
        std::vector<double> mst_weights;
        mst_weights.reserve(mst.size());
        for (auto [u,v] : mst) {
            mst_weights.push_back(m_graph.getEdgeDistance(u,v));
        }

        ultrametric.buildFromMST(m_graph.getNumNodes(), mst, mst_weights);

        assert(ultrametric.root != -1);

    }


    void RaeckeCKR::setGraph(const Graph& g) {
        m_graph = g;
    }

    const Graph& RaeckeCKR::getGraph() const {
        return m_graph;
    }



    TreeNode* RaeckeCKR::getTree(Graph& g) {
        const int n = g.getNumNodes();
        if (n == 0) return nullptr;


        std::mt19937 rng(std::random_device{}());
        m_levels.clear();

        // --- (0) preprocess ultrametric from MST ---
        preprocess();

        // ---- (1) prepare node pointers for the finest level ----
        std::vector<TreeNode*> prev_nodes(n);
        for (int v = 0; v < n; ++v) {
            prev_nodes[v] = new TreeNode(v);
            prev_nodes[v]->members = {v};  // keep original vertex IDs here
        }

        // ---- (2) choose logarithmic set of Δ-scales ----
        std::vector<double> scales;
        double Delta = diameter;
        for (; Delta >= diameter / (2.0 * n); Delta /= 8.0)
            scales.push_back(Delta);

        if (debug)
            std::cout << "[MendelScaling] Using " << scales.size()
                  << " quotient scales (log n levels)\n";

        // ---- (3) build hierarchical levels ----
        TreeNode* root = nullptr;
        std::vector<int> current_centers(n);
        std::iota(current_centers.begin(), current_centers.end(), 0);

        for (double Delta : scales) {
            QuotientLevel Q = build_quotient_graph_with_map(g, ultrametric, Delta);
            if (Q.Gq.getNumNodes() <= 1) {
                // if single cluster: make that the root once at the end
                continue;
            }

            m_levels.emplace_back();
            CKRLevel& L = m_levels.back();
            build_ckr_level(Q.Gq, Delta, L, rng);

            // parent nodes at this Δ-level
            std::unordered_map<int, TreeNode*> center_to_parent;
            center_to_parent.reserve(L.centers.size()*2);

            // --- group vertices by owner into new TreeNodes (as in your original code) ---
            std::unordered_map<int, TreeNode*> center_to_node;
            // Efficiently map each prev_node (cluster from previous finer step) to a quotient vertex id
            // We use the first original member of that cluster as representative:
            for (TreeNode* childCluster : prev_nodes) {
                if (!childCluster) continue;
                int rep = childCluster->members.empty() ? -1 : childCluster->members[0];
                if (rep < 0) continue;

                // find qid for representative and then its owner center c
                int qid = Q.sigma_compact_of_v[rep];
                int c = L.owner[qid];
                if (c == -1) continue;  // should be rare; skip if not assigned (due to R cutoff)

                // parent TreeNode for this center
                TreeNode* parent = nullptr;
                auto it = center_to_parent.find(c);
                if (it == center_to_parent.end()) {
                    parent = new TreeNode(L.centers[c]);
                    parent->radius = Delta;
                    center_to_parent[c] = parent;
                } else parent = it->second;

                // Attach childCluster under parent; merge members
                parent->children.push_back(childCluster);
                childCluster->parent = parent;
                parent->members.insert(parent->members.end(),
                                       childCluster->members.begin(),
                                       childCluster->members.end());
            }

            // Prepare for next iteration (finer scale): current parents become the "prev_nodes"
            prev_nodes.clear();
            prev_nodes.reserve(center_to_parent.size());
            for (auto& [_, node] : center_to_parent)
                prev_nodes.push_back(node);
        }

        // Final root: if multiple parents remain, wrap them under one root; else take the only one.
        if (prev_nodes.empty()) {
            // Degenerate case: all collapsed early → build a single root that contains all vertices
            root = new TreeNode(0);
            root->members.resize(n);
            std::iota(root->members.begin(), root->members.end(), 0);
        } else if (prev_nodes.size() == 1) {
            root = prev_nodes.front();
        } else {
            root = new TreeNode(-1);
            root->radius = (scales.empty()? 0.0 : scales.front());
            for (auto* ch : prev_nodes) {
                root->children.push_back(ch);
                ch->parent = root;
                root->members.insert(root->members.end(), ch->members.begin(), ch->members.end());
            }
        }

        // sanity: all original vertices must be present exactly once
        // (you can add a debug assert that counts coverage here)
        return root;
    }




    void RaeckeCKR::computeRLoads(TreeNode *t, Graph &g) {
    // iterate through the tree and store the rloads into the adjacency list
    std::queue<TreeNode*> q;
        q.push(t);

        while (!q.empty()) {
            TreeNode* node = q.front();
            q.pop();

            // --- 3️⃣ Process each child: represents a cut S_child | V\S_child ---
            for (TreeNode* child : node->children) {
                // Add child to traversal queue
                q.push(child);

                const std::vector<int>& clusterVertices = child->members;
                if (clusterVertices.empty()) continue;

                // Build set for fast lookup
                std::vector<char> S(m_graph.getNumNodes(), 0);
                for (int v : clusterVertices) S[v] = 1;


                // --- 4️⃣ Compute total cut capacity of this child cluster ---
                double cut = 0.0;
                for (int u : clusterVertices) {
                    for (int v : m_graph.neighbors(u)) {
                        if (!S[v])  // boundary edge
                            cut += m_graph.getEdgeCapacity(u, v);
                    }
                }
                if (cut <= 1e-12) cut = 1e-12;  // avoid zero-division

                // --- 5️⃣ Choose representative vertices (simple heuristic) ---
                int repParent = node->members.empty() ? clusterVertices[0] : node->members[0];
                int repChild  = clusterVertices[0];

                // --- 6️⃣ Compute shortest path between parent and child reps ---
                auto path = m_graph.getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                // --- 7️⃣ Update edge r-loads along that path ---
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int u = path[i];
                    int v = path[i + 1];
                    double cap = m_graph.getEdgeCapacity(u, v);
                    if (cap <= 1e-12) cap = 1e-12;

                    double delta = cut / cap;
                    addLoadToEdge(u, v, delta);

                }
            }
        }
    }



    // TODO: this is for now a naive not optimized, for each newly added load => O(n) linear possible -> could be reduced to O(log) at least with binary search tree
    void RaeckeCKR::addLoadToEdge(int u, int v, double load) {
        edge2Load[{u, v}] += load;
        edge2Load[{v, u}] += load;
    }


    double RaeckeCKR::getMaxRload() {
        double max_load = 0.0;

        for (const auto& [edge, load] : edge2Load) {
            if (load > max_load) {
                max_load = load;
            }
        }

        return max_load;
    }

    void RaeckeCKR::computeNewDistances(Graph &g) {
        constexpr double EPS = 1e-12;

        // 1) Build total_r(e) for all edges, default 0
        std::unordered_map<std::pair<int,int>, double> total_r;
        total_r.reserve(g.getNumEdges());

        for (int u = 0; u < g.getNumNodes(); ++u) {
            for (int v : g.neighbors(u)) {
                int a = std::min(u,v), b = std::max(u,v);
                // only one direction
                if (u != a) continue;
                auto it = edge2Load.find({a,b});
                total_r[{a,b}] = (it == edge2Load.end()) ? 0.0 : it->second;
            }
        }

        // 2) softmax with shift
        double max_r = 0.0;
        for (auto& [e,r] : total_r) max_r = std::max(max_r, r);

        double sumExp = 0.0;
        for (auto& [e,r] : total_r) sumExp += std::exp(r - max_r);
        if (sumExp <= 0.0 || !std::isfinite(sumExp)) sumExp = 1.0;

        // 3) distances
        double min_d = std::numeric_limits<double>::infinity();
        std::unordered_map<std::pair<int,int>, double> newDist;
        newDist.reserve(total_r.size());

        for (auto& [e,r] : total_r) {
            auto [a,b] = e;
            double cap = g.getEdgeCapacity(a,b);
            if (cap < EPS) cap = EPS;
            double d = (std::exp(r - max_r) / cap) / sumExp;
            newDist[e] = d;
            if (d < min_d) min_d = d;
        }
        if (min_d < EPS) min_d = EPS;

        for (auto& [e,d] : newDist) {
            auto [a,b] = e;
            double norm = d / min_d;
            if (norm < 1.0) norm = 1.0;
            g.updateEdgeDistance(a,b, norm);
            g.updateEdgeDistance(b,a, norm); // mirror, if your graph stores both arcs
        }
    }



    //// Flow Storage methods ////
    void FlowStorage::init(const Graph& g) {
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

        void FlowStorage::addFlow(int u, int v, int s, double flow) {
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

                addFlow(lo, s, flow);

        }

    // -----------------------------------------------------------------------------
    // ⚡ Optimized addFlow (manual binary search, cache-friendly)
    // -----------------------------------------------------------------------------
    #if defined(__GNUC__) || defined(__clang__)
        __attribute__((always_inline))
        #elif defined(_MSC_VER)
        __forceinline
        #endif
        inline void FlowStorage::addFlow(int edge_id, int u, double flow) noexcept {
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

        void FlowStorage::scaleDownFlow() {
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


    ///// CKRTransform methods ////
    void CKRTransform::init(const Graph& g) {
        // fix a node for the storing the lienar oblivious routing
        this->flow_storage.init(g);
    }

    inline void CKRTransform::addFlow(int u, int v, int s, int t, double flow) {
        flow_storage.addFlow(u, v, s, flow);
    }

    void CKRTransform::addTree(const TreeNode* root, double lambda, Graph& graph) {
        if (!root) return;
        std::queue<const TreeNode*> q;
        q.push(root);

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


}