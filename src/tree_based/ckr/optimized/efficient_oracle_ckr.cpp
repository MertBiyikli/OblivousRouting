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



/*
 * Preprocess the graph to build the ultrametric tree from MST.
 * The idea is to first compute an MST (using Kruskal) for computing an minimum spanning tree.
 * Then, we build the ultrametric tree from the MST edges and their weights.
 * The goal is to preprocess the algorithm such that the initial solution(MST) n-1 approximates the optimal solution.
 * This will then be used to construct the ultra tree metric. See paper: Har-Peled, Sariel, and Manor Mendel. "Fast construction of nets in low dimensional metrics, and their applications."
 */
void EfficientCKR::preprocess() {
    // -- compute MST using Kruskals Algo
    RandomMST mst_oracle;
    mst_oracle.setGraph(*g_ptr);
    auto mst = mst_oracle.build_mst();


    // ----  Collect MST edge weights ----
    std::vector<double> mst_weights;
    mst_weights.reserve(mst.size());
    for (auto [u,v] : mst) {
        mst_weights.push_back(g_ptr->getEdgeDistance(u,v));
    }

    ultrametric.buildFromMST(g_ptr->getNumNodes(), mst, mst_weights);

    assert(ultrametric.root != -1);
}


/*
 * Build the CKR tree from the preprocessed ultrametric tree and the given graph.
 * The idea is to use the ultrametric tree to guide the partitioning of the graph at different scales (Δ).
 * For each scale, we build a quotient graph and then apply the CKR partitioning algorithm to create clusters.
 * These clusters are then organized into a hierarchical tree structure.
 *
 * There are two major efficiency improvements here:
 * 1) We use the ultrametric tree efficiently build the subgraphs that are partitioned on each level.
 * 2) Instead of naive O(diameter) scales, we use logarithmic number of scales based on Mendel and Schwob's theoretical results.
 */
std::shared_ptr<TreeNode> EfficientCKR::getTree() {
    if (g_ptr == nullptr)
        throw std::runtime_error("CKR: graph pointer is null");

    const int n = g_ptr->getNumNodes();
    if (n == 0) return nullptr;

    m_levels.clear();
    // --- (0) preprocess ultrametric from MST ---
    preprocess();


    // ---- (1) prepare node pointers for the finest level ----
    std::vector<std::shared_ptr<TreeNode>> prev_nodes(n);
    for (int v = 0; v < n; ++v) {
        prev_nodes[v] = std::make_shared<TreeNode>();
        prev_nodes[v]->id = v;
        prev_nodes[v]->members = {v};  // keep original vertex IDs here
    }

    // ---- (2) choose logarithmic set of Δ-scales ----
    computeScales();

    // ---- (3) build hierarchical levels ----
    std::shared_ptr<TreeNode> root;

    MendelScaling::QuotientConstruction<Graph_csr> qc;
    qc.preprocessEdges(dynamic_cast<const Graph_csr &>(*g_ptr));


    // ---- (4) for each scale, build the CKR level and link to previous level ----
    for (double Delta : scales) {
        MendelScaling::QuotientLevel<Graph_csr> Q;

        Q = qc.constructQuotientGraph(ultrametric, Delta,dynamic_cast<Graph_csr &>(*g_ptr));
        if (Q.Gq.getNumNodes() <= 1) {
            // if single cluster: make that the root once at the end
            continue;
        }

        m_levels.emplace_back();
        CKRLevel& L = m_levels.back();
        computeLevelPartition(Q, Delta, L);

        buildTreeLevel(prev_nodes, Q, L, Delta);

    }

    finishTree(root, prev_nodes);
    // sanity: all original vertices must be present exactly once
    // (you can add a debug assert that counts coverage here)
    return root;
}


void EfficientCKR::computeScales() {
    const int n = g_ptr->getNumNodes();
    scales.clear();
    if (n == 0) return;
    double min_weight = std::numeric_limits<double>::max();
    for (int e = 0; e<g_ptr->getNumEdges(); e++) {
        double w = g_ptr->getEdgeDistance(e);
        if (w < min_weight) {
            min_weight = w;
        }
    }

    // compute the new diameter after preprocessing
    if (diameter == 0) {
        throw std::runtime_error("Graph has zero diameter; cannot build CKR tree.");
    }
    double d = diameter * 2 * n;
    for (; d >= diameter / (2.0 * n) || d > min_weight; d /= 8.0)
        scales.push_back(d);

    // sort the scales in increasing order
    std::reverse(scales.begin(), scales.end());
}


void EfficientCKR::computeLevelPartition(MendelScaling::QuotientLevel<Graph_csr> &Q, double Delta, CKRLevel &L) {

    auto start_time = std::chrono::high_resolution_clock::now();
    build_ckr_level(Q.Gq, Delta, L);
    auto end_time = std::chrono::high_resolution_clock::now();
    this->pure_oracle_running_times.emplace_back(std::chrono::duration<double, std::milli>(end_time - start_time).count());

}


void EfficientCKR::buildTreeLevel( std::vector<std::shared_ptr<TreeNode> > &prev_nodes, const MendelScaling::QuotientLevel<Graph_csr> &Q, const CKRLevel &L, const double Delta) {
    // parent nodes at this Δ-level
    std::unordered_map<int, std::shared_ptr<TreeNode>> center_to_parent;
    center_to_parent.reserve(L.centers.size()*2);


    // Efficiently map each prev_node (cluster from previous finer step) to a quotient vertex id
    // We use the first original member of that cluster as representative:
    for (const auto& childCluster : prev_nodes) {
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



void EfficientCKR::finishTree(std::shared_ptr<TreeNode>& root, const std::vector<std::shared_ptr<TreeNode> > &prev_nodes) {
    if (!g_ptr) {
        return;
    }
    // Final root: if multiple parents remain, wrap them under one root; else take the only one.
    if (prev_nodes.empty()) {
        // Degenerate case: all collapsed early → build a single root that contains all vertices
        root = std::make_shared< TreeNode>();
        root->id = 0;
        root->members.resize(g_ptr->getNumNodes());
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
}




std::vector<int> EfficientCKR::build_ckr_level(const Graph_csr &g, double Delta, CKRLevel &L) {
    CKRPartition ckr;

    const int n = g.getNumNodes();
    std::vector<int> X(n);
    std::iota(X.begin(), X.end(), 0);
    auto clusters = ckr.computePartition(g, X, Delta, L);

    return clusters;
}


/*

void EfficientCKR::computeRLoads(std::shared_ptr<TreeNode> t) {
        assert(t != nullptr);
        assert(g_ptr != nullptr);
        // clear current rloads

    std::fill(rload_current.begin(), rload_current.end(), 0.0);
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

                    int e = g_ptr->getEdgeId(u, v);
                    double rload = rload_current[e];
                    double cap = g_ptr->getEdgeCapacity(u, v);
                    if (cap <= 1e-12) cap = 1e-12;

                    rload += cut / cap;
                    rload_current[e] = rload;

                    int anti_e = g_ptr->getEdgeId(v, u);
                    rload_current[anti_e] = rload;

                }
            }
        }
    }




double EfficientCKR::getMaxRload() const {
    assert(rload_current.size() > 0);

    double max_load = 0.0;
    for (const auto& load : rload_current) {
        if (load > max_load) {
            max_load = load;
        }
    }
    return max_load;
}


void EfficientCKR::computeNewDistances(double lambda) {
    assert(g_ptr != nullptr);
    constexpr double EPS = 1e-12;

    addCurrentLoad(lambda);

    double max_r = 0.0;
    for (auto& r : rload_total) max_r = std::max(max_r, r);

    double sumExp = 0.0;
    for (auto& r : rload_total) sumExp += std::exp(r - max_r);
    if (sumExp <= 0.0 || !std::isfinite(sumExp)) sumExp = 1.0;

    // 3) distances
    double min_d = std::numeric_limits<double>::infinity();
    std::vector<double> newDist;
    newDist.reserve(rload_total.size());

    for (int e = 0; e < g_ptr->getNumEdges(); ++e) {
        double r = rload_total[e];

        double cap = g_ptr->getEdgeCapacity(e);
        if (cap < EPS) cap = EPS;
        double d = (std::exp(r - max_r) / cap) / sumExp;
        newDist[e] = d;
        if (d < min_d) min_d = d;
    }
    if (min_d < EPS) min_d = EPS;

    for (int e = 0; e < g_ptr->getNumEdges(); ++e) {
        double d = newDist[e];
        double norm = d / min_d;
        if (norm < 1.0) norm = 1.0;
        g_ptr->updateEdgeDistance(e, norm);
        //g.updateEdgeDistance(b,a, norm); // mirror, if your graph stores both arcs
    }

}

void EfficientCKR::addCurrentLoad(double lambda) {
    assert(g_ptr != nullptr);
    for (int e = 0; e < g_ptr->getNumEdges(); ++e) {
        rload_total[e] += rload_current[e] * lambda;
    }
}
*/