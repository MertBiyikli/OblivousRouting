//
// Created by Mert Biyikli on 24.11.25.
//

#include "efficient_oracle_ckr.h"
#include "../utils/quotient_graph.h"
#include "../utils/ultrametric_tree.h"
#include "../../random_mst/mst.h"
#include <cassert>
#include <random>
#include "../../raecke_tree.h"
#include "../../../datastructures/GraphCSR.h"



/*
 * Preprocess the graph to build the ultrametric tree from MST.
 * The idea is to first compute an MST (using Kruskal) for computing an minimum spanning tree.
 * Then, we build the ultrametric tree from the MST edges and their weights.
 * The goal is to preprocess the algorithm such that the initial solution(MST) n-1 approximates the optimal solution.
 * This will then be used to construct the ultra tree metric. See paper: Har-Peled, Sariel, and Manor Mendel. "Fast construction of nets in low dimensional metrics, and their applications."
 */
void EfficientCKR::preprocess() {
    // -- compute MST using Kruskals Algo
    RandomMST mst_oracle(g);
    //mst_oracle.setGraph(g);
    auto mst = mst_oracle.build_mst();


    // ----  Collect MST edge weights ----
    std::vector<double> mst_weights;
    mst_weights.reserve(mst.size());
    for (auto [u,v] : mst) {
        mst_weights.push_back(g.getEdgeDistance(u,v));
    }

    ultrametric.buildFromMST(g.getNumNodes(), mst, mst_weights);

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

    const int n = g.getNumNodes();
    if (n == 0) return nullptr;

    diameter = g.GetDiameter();

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

    MendelScaling::QuotientConstruction qc;
    qc.preprocessEdges(g);


    // ---- (4) for each scale, build the CKR level and link to previous level ----
    for (double Delta : scales) {
        MendelScaling::QuotientLevel Q;

        Q = qc.constructQuotientGraph(ultrametric, Delta,g);
        if (Q.Gq->getNumNodes() <= 1) {
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
    const int n = g.getNumNodes();
    scales.clear();
    if (n == 0) return;
    double min_weight = std::numeric_limits<double>::max();
    for (int e = 0; e< g.getNumEdges(); e++) {
        double w = g.getEdgeDistance(e);
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


void EfficientCKR::computeLevelPartition(MendelScaling::QuotientLevel &Q, double Delta, CKRLevel &L) {
    auto start_time = std::chrono::high_resolution_clock::now();
    build_ckr_level(*Q.Gq, Delta, L);
    auto end_time = std::chrono::high_resolution_clock::now();
    //this->pure_oracle_running_times.emplace_back(std::chrono::duration<double, std::milli>(end_time - start_time).count());

}


void EfficientCKR::buildTreeLevel( std::vector<std::shared_ptr<TreeNode> > &prev_nodes, const MendelScaling::QuotientLevel &Q, const CKRLevel &L, const double Delta) {
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

    // Final root: if multiple parents remain, wrap them under one root; else take the only one.
    if (prev_nodes.empty()) {
        // Degenerate case: all collapsed early → build a single root that contains all vertices
        root = std::make_shared< TreeNode>();
        root->id = 0;
        root->members.resize(g.getNumNodes());
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




std::vector<int> EfficientCKR::build_ckr_level(const IGraph &g, double Delta, CKRLevel &L) {
    CKRPartition ckr;

    const int n = g.getNumNodes();
    std::vector<int> X(n);
    std::iota(X.begin(), X.end(), 0);
    auto clusters = ckr.computePartition(g, X, Delta, L);

    return clusters;
}


void EfficientCKR::updateEdgeDistances(const std::vector<double> &distances) {
    for (int e = 0; e < g.getNumEdges(); ++e) {
        g.updateEdgeDistance(e, distances[e]);
    }
}
