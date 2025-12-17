//
// Created by Mert Biyikli on 24.11.25.
//

#include "efficient_oracle_ckr.h"
#include "../utils/quotient_graph.h"
#include "../utils/ultrametric_tree.h"
#include "../../random_mst/mst.h"
#include <cassert>
#include <random>
#include <google/protobuf/arena_cleanup.h>

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
        prev_nodes[v]->center = v;
    }

    // ---- (2) choose logarithmic set of Δ-scales ----
    computeScales();

    // ---- (3) build hierarchical levels ----
    std::shared_ptr<TreeNode> root = std::make_shared<TreeNode>();

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
    CKRPartition ckr;
    build_ckr_level(ckr, *Q.Gq, Delta, L);
    auto qid = ckr.build_qid_to_rep(Q, g.getNumNodes(), Q.Gq->getNumNodes());
    ckr.build_cluster_of_qid(Q, qid, L);

    auto end_time = std::chrono::high_resolution_clock::now();
    //this->pure_oracle_running_times.emplace_back(std::chrono::duration<double, std::milli>(end_time - start_time).count());

}


std::shared_ptr<TreeNode>
EfficientCKR::buildTreeLevel(
std::vector<std::shared_ptr<TreeNode>>& prev_nodes,
        const MendelScaling::QuotientLevel& Q,
        const CKRLevel& L,
        const double Delta
) {
    // ------------------------------------------------------------------
    // Step 0: Map quotient nodes → representative original vertex
    // ------------------------------------------------------------------
    const int nq = Q.Gq->getNumNodes();
    std::vector<int> qid_to_rep(nq, -1);

    for (int v = 0; v < g.getNumNodes(); ++v) {
        int qid = Q.sigma_compact_of_v[v];
        if (qid >= 0 && qid < nq && qid_to_rep[qid] == -1)
            qid_to_rep[qid] = v;
    }

    // ------------------------------------------------------------------
    // Step 1: Create parent nodes (one per CKR cluster)
    // ------------------------------------------------------------------
    std::vector<std::shared_ptr<TreeNode>> parents(L.centers.size());

    for (size_t c = 0; c < L.centers.size(); ++c) {
        auto parent = std::make_shared<TreeNode>();

        const int center_qid = L.centers[c];
        int center_rep = -1;

        if (center_qid >= 0 && center_qid < nq)
            center_rep = qid_to_rep[center_qid];

        // Fallback (should almost never happen)
        if (center_rep == -1 && !prev_nodes.empty())
            center_rep = prev_nodes.front()->center;

        parent->center = center_rep;
        parent->parent.reset();
        parent->children.clear();
        parent->members.clear();

        parents[c] = parent;
    }

    // ------------------------------------------------------------------
    // Step 2: Attach lower-level nodes to parents
    //         (according to CKR assignment)
    // ------------------------------------------------------------------
    for (const auto& child : prev_nodes) {
        if (!child) continue;

        int child_center = child->center;
        int qid = Q.sigma_compact_of_v[child_center];
        int cluster = L.cluster_of_qid[qid];

        auto& parent = parents[cluster];

        child->parent = parent;
        parent->children.push_back(child);

        // merge members upward
        parent->members.insert(
            parent->members.end(),
            child->members.begin(),
            child->members.end()
        );
    }

    // ------------------------------------------------------------------
    // Step 3: Deduplicate members (important!)
    // ------------------------------------------------------------------
    for (auto& p : parents) {
        auto& M = p->members;
        std::sort(M.begin(), M.end());
        M.erase(std::unique(M.begin(), M.end()), M.end());
    }

    // ------------------------------------------------------------------
    // Step 4: If this is the topmost level, return virtual root
    // ------------------------------------------------------------------
    if (L.R == scales.front()) {
        auto root = std::make_shared<TreeNode>();
        root->parent.reset();
        root->children.clear();
        root->members.clear();

        // Choose a stable center
        root->center = parents.front()->center;

        for (auto& p : parents) {
            p->parent = root;
            root->children.push_back(p);
            root->members.insert(
                root->members.end(),
                p->members.begin(),
                p->members.end()
            );
        }

        std::sort(root->members.begin(), root->members.end());
        root->members.erase(
            std::unique(root->members.begin(), root->members.end()),
            root->members.end()
        );

        return root;
    }

    // Otherwise return a dummy aggregator node (used by caller)
    auto aggregator = std::make_shared<TreeNode>();
    aggregator->parent.reset();
    aggregator->children = parents;
    aggregator->members.clear();

    aggregator->center = parents.front()->center;

    for (auto& p : parents) {
        p->parent = aggregator;
        aggregator->members.insert(
            aggregator->members.end(),
            p->members.begin(),
            p->members.end()
        );
    }

    std::sort(aggregator->members.begin(), aggregator->members.end());
    aggregator->members.erase(
        std::unique(aggregator->members.begin(), aggregator->members.end()),
        aggregator->members.end()
    );

    return aggregator;
}



void EfficientCKR::finishTree(std::shared_ptr<TreeNode>& root, const std::vector<std::shared_ptr<TreeNode> > &prev_nodes) {

    // Final root: if multiple parents remain, wrap them under one root; else take the only one.
    if (prev_nodes.empty()) {
        // Degenerate case: all collapsed early → build a single root that contains all vertices
        root = std::make_shared< TreeNode>();
        root->id = 0;
        root->members.resize(g.getNumNodes());
        std::iota(root->members.begin(), root->members.end(), 0);
        root->center = (!root->members.empty() ? root->members[0] : -1 );
    } else if (prev_nodes.size() == 1) {
        root = prev_nodes.front();
    } else {
        root =std::make_shared< TreeNode>();
        //root->id = -1;
        root->radius = (scales.empty()? 0.0 : scales.front());
        root->center = prev_nodes.front()->center;
        for (auto ch : prev_nodes) {
            root->children.push_back(ch);
            ch->parent = root;
            root->members.insert(root->members.end(), ch->members.begin(), ch->members.end());
        }
    }

    cleanUpTree(root);
}


void EfficientCKR::cleanUpTree(std::shared_ptr<TreeNode>& node) {
    if (!node) return;

    // first clean children
    for (auto& ch : node->children) cleanUpTree(ch);

    // then rebuild children list safely
    std::vector<std::shared_ptr<TreeNode>> new_children;
    new_children.reserve(node->children.size());

    for (auto& child : node->children) {
        if (!child) continue;

        if (same_members(node->members, child->members)) {
            // absorb child: adopt grandchildren
            for (auto& gc : child->children) {
                if (!gc) continue;
                gc->parent = node;
                new_children.push_back(gc);
            }
            child->children.clear();
            // child gets dropped
        } else {
            new_children.push_back(child);
        }
    }

    node->children.swap(new_children);
}




std::vector<int> EfficientCKR::build_ckr_level( CKRPartition& ckr, const IGraph &g, double Delta, CKRLevel &L) {

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
