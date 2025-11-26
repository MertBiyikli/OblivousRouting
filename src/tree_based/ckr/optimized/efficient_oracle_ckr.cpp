//
// Created by Mert Biyikli on 24.11.25.
//

#include "efficient_oracle_ckr.h"
#include "../utils/quotient_graph.h"
#include "../utils/ultrametric_tree.h"
#include "../../random_mst/mst.h"
#include <cassert>
#include <random>

#include "../../../datastructures/graph_csr.h"



void EfficientCKR::preprocess() {
        // -- compute MST using Kruskals Algo
        RandomMST mst_oracle;
        mst_oracle.setGraph(*g_ptr);
        std::vector<std::pair<int, int> > mst = mst_oracle.build_mst();


        // ----  Collect MST edge weights ----
        std::vector<double> mst_weights;
        mst_weights.reserve(mst.size());
        for (auto [u,v] : mst) {
            mst_weights.push_back(g_ptr->getEdgeDistance(u,v));
        }

        ultrametric.buildFromMST(g_ptr->getNumNodes(), mst, mst_weights);

        assert(ultrametric.root != -1);

    }




    std::shared_ptr<TreeNode> EfficientCKR::getTree() {
        const int n = g_ptr->getNumNodes();
        if (n == 0) return nullptr;


        std::mt19937 rng(std::random_device{}());
        m_levels.clear();

        // --- (0) preprocess ultrametric from MST ---
        preprocess();

        if (debug ) {
            // print ultrametric tree
            std::cout << "[MendelScaling] Ultrametric Tree: " << std::endl;
            for (int i = 0; i < ultrametric.T.size(); i++) {
                std::cout << "Node " << i << " | Gamma: " << ultrametric.T[i].Gamma << " | Children: ";
                for (auto& c : ultrametric.T[i].child) {
                    std::cout << c << " ";
                }
                std::cout << std::endl;
            }
        }

        // ---- (1) prepare node pointers for the finest level ----
        std::vector<std::shared_ptr<TreeNode>> prev_nodes(n);
        for (int v = 0; v < n; ++v) {
            prev_nodes[v] = std::make_shared<TreeNode>();
            prev_nodes[v]->id = v;
            prev_nodes[v]->members = {v};  // keep original vertex IDs here
        }

        // ---- (2) choose logarithmic set of Δ-scales ----
        // but first identify minimum edge weight
        double min_weight = std::numeric_limits<double>::max();
        for (int e = 0; e<g_ptr->getNumEdges(); e++) {
            double w = g_ptr->getEdgeDistance(e);
            if (w < min_weight) {
                min_weight = w;
            }
        }

        std::vector<double> scales;
        // compute the new diameter after preprocessing

        if (diameter == 0) {
         throw std::runtime_error("Graph has zero diameter; cannot build CKR tree.");
        }
        double d = diameter * 2 * n;
        for (; d >= diameter / (2.0 * n) || d > min_weight; d /= 8.0)
            scales.push_back(d);

        // sort the scales in increasing order
        std::reverse(scales.begin(), scales.end());

        if (debug)
            std::cout << "[MendelScaling] Using " << scales.size()
                  << " quotient scales (log n levels)\n";

        // ---- (3) build hierarchical levels ----
        std::shared_ptr<TreeNode> root;
        std::vector<int> current_centers(n);
        std::iota(current_centers.begin(), current_centers.end(), 0);

        for (double Delta : scales) {
            if (debug )
                std::cout << "[MendelScaling] Building level for Δ = " << Delta << "\n";


            MendelScaling::QuotientLevel<Graph_csr> Q = build_quotient_graph_with_map(*g_ptr, ultrametric, Delta);
            if (Q.Gq.getNumNodes() <= 1) {
                // if single cluster: make that the root once at the end
                continue;
            }

            m_levels.emplace_back();
            CKRLevel& L = m_levels.back();
            auto start_time = std::chrono::high_resolution_clock::now();
            build_ckr_level(Q.Gq, Delta, L);
            auto end_time = std::chrono::high_resolution_clock::now();
            this->pure_oracle_running_times.emplace_back(std::chrono::duration<double, std::milli>(end_time - start_time).count());


            // for debugging purpose print all information,e .g. Quotient graph and the computed partitons
            if (debug ) {
                std::cout << "[MendelScaling] Δ = " << Delta
                          << ", Quotient graph nodes: " << Q.Gq.getNumNodes()
                          << ", edges: " << Q.Gq.getNumEdges() << "\n";
                std::cout << "[MendelScaling] CKR centers at this level: ";
                for (auto c : L.centers)
                    std::cout << c << " ";
                std::cout << "\n";

                std::cout << "With labels: " << std::endl;
                for (int v = 0; v < Q.Gq.getNumNodes(); ++v) {
                    std::cout << "Vertex " << v << " assigned to center " << L.owner[v] << std::endl;
                }

                std::cout << "Current Quotient Graph : " << std::endl;
                Q.Gq.print();
            }

            // parent nodes at this Δ-level
            std::unordered_map<int, std::shared_ptr<TreeNode>> center_to_parent;
            center_to_parent.reserve(L.centers.size()*2);

            // --- group vertices by owner into new TreeNodes (as in your original code) ---
            std::unordered_map<int, std::shared_ptr<TreeNode>> center_to_node;
            // Efficiently map each prev_node (cluster from previous finer step) to a quotient vertex id
            // We use the first original member of that cluster as representative:
            for (auto childCluster : prev_nodes) {
                if (!childCluster) continue;
                int rep = childCluster->members.empty() ? -1 : childCluster->members[0];
                if (rep < 0) continue;

                // find qid for representative and then its owner center c
                int qid = Q.sigma_compact_of_v[rep];
                int c = L.owner[qid];
                if (c == -1) continue;  // should be rare; skip if not assigned (due to R cutoff)

                // parent TreeNode for this center
                std::shared_ptr<TreeNode> parent;
                auto it = center_to_parent.find(c);
                if (it == center_to_parent.end()) {
                    parent = std::make_shared<TreeNode>();
                    parent->id = L.centers[c];
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
            root = std::make_shared< TreeNode>();
            root->id = 0;
            root->members.resize(n);
            std::iota(root->members.begin(), root->members.end(), 0);
        } else if (prev_nodes.size() == 1) {
            root = prev_nodes.front();
        } else {
            root =std::make_shared< TreeNode>();
            root->id = -1;
            root->radius = (scales.empty()? 0.0 : scales.front());
            for (auto ch : prev_nodes) {
                root->children.push_back(ch);
                ch->parent = root;
                root->members.insert(root->members.end(), ch->members.begin(), ch->members.end());
            }
        }

        // sanity: all original vertices must be present exactly once
        // (you can add a debug assert that counts coverage here)
        return root;
    }


    std::vector<int> EfficientCKR::build_ckr_level(const Graph_csr &g, double Delta, CKRLevel &L) {
        CKRPartition ckr;

        const int n = g.getNumNodes();
        std::vector<int> X(n);
        std::iota(X.begin(), X.end(), 0);
        auto clusters = ckr.computePartition(g, X, Delta, L);

        return clusters;
    }




    void EfficientCKR::computeRLoads(std::shared_ptr<TreeNode> t, int tree_index) {
        // --- 1️⃣ Ensure enough space for r-load maps ---
        if (edge2Load.size() <= tree_index)
            edge2Load.resize(tree_index + 1);

        auto& local_edge2Load = edge2Load[tree_index];
        local_edge2Load.clear();

        // --- 2️⃣ BFS or DFS traversal (we'll use queue for clarity) ---
    // iterate through the tree and store the rloads into the adjacency list
    std::queue<std::shared_ptr<TreeNode>> q;
        q.push(t);

        while (!q.empty()) {
            std::shared_ptr<TreeNode> node = q.front();
            q.pop();

            // --- 3️⃣ Process each child: represents a cut S_child | V\S_child ---
            for (auto child : node->children) {
                // Add child to traversal queue
                q.push(child);


                // if node has the same nodes as child, skip
                if (node->members.size() == child->members.size()) {
                    continue;
                }
                const std::vector<int>& clusterVertices = child->members;
                if (clusterVertices.empty()) continue;

                // Build set for fast lookup
                std::vector<char> S(g_ptr->getNumNodes(), 0);
                for (int v : clusterVertices) S[v] = 1;


                // --- 4️⃣ Compute total cut capacity of this child cluster ---
                double cut = 0.0;
                for (int u : clusterVertices) {
                    for (auto&  v : g_ptr->neighbors(u)) {
                        if (!S[v])  // boundary edge
                            cut += g_ptr->getEdgeCapacity(u, v);
                    }
                }
                if (cut <= 1e-12) cut = 1e-12;  // avoid zero-division

                // --- 5️⃣ Choose representative vertices (simple heuristic) ---
                int repParent = node->members.empty() ? clusterVertices[0] : node->members[0];
                int repChild  = clusterVertices[0];

                // --- 6️⃣ Compute shortest path between parent and child reps ---
                auto path = g_ptr->getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                // --- 7️⃣ Update edge r-loads along that path ---
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int u = path[i];
                    int v = path[i + 1];

                    double rload = local_edge2Load[{u, v}];
                    double cap = g_ptr->getEdgeCapacity(u, v);
                    if (cap <= 1e-12) cap = 1e-12;

                    rload += cut / cap;
                    local_edge2Load[{u, v}] = rload;
                    local_edge2Load[{v, u}] = rload;

                }
            }
        }
    }




    double EfficientCKR::getMaxRload(int tree_index) const {
        double max_load = 0.0;

        for (const auto& [edge, load] : edge2Load[tree_index]) {
            if (load > max_load) {
                max_load = load;
            }
        }

        return max_load;
    }

    void EfficientCKR::computeNewDistances() {
        constexpr double EPS = 1e-12;

        // 1) Build total_r(e) for all edges, default 0
        std::unordered_map<std::pair<int,int>, double> total_r;
        total_r.reserve(g_ptr->getNumEdges());

        double max_total_r = -std::numeric_limits<double>::infinity();

        for (int u = 0; u < g_ptr->getNumNodes(); ++u) {
            for (auto& v : g_ptr->neighbors(u)) {
                for (int tree_index = 0; tree_index < edge2Load.size(); ++tree_index) {
                    auto& edge2Load_map = edge2Load[tree_index];
                    auto it = edge2Load_map.find({u,v});
                    if (it != edge2Load_map.end()) {
                        total_r[{u,v}] += it->second * m_lambdas[tree_index];
                    }
                }
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
            double cap = g_ptr->getEdgeCapacity(a,b);
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
            g_ptr->updateEdgeDistance(a,b, norm);
            //g.updateEdgeDistance(b,a, norm); // mirror, if your graph stores both arcs
        }
    }


