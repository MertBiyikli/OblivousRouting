//
// Created by Mert Biyikli on 24.11.25.
//

#ifndef OBLIVIOUSROUTING_EFFICIENT_CKR_H
#define OBLIVIOUSROUTING_EFFICIENT_CKR_H

#include "../../../solver/solver.h"
#include "../../../datastructures/GraphCSR.h"
#include "../ckr_tree_decomposer.h"
#include "../raecke_ckr_transform.h"
#include "../../raecke_oracle_iteration.h"
#include "../utils/ultrametric_tree.h"
#include "../ckr_tree_decomposer.h"
#include "../utils/quotient_graph.h"
#include "../../raecke_tree.h"
#include "../../efficient_raecke_base.h"
#include <chrono>

/*
    * Optimized CKR implementation with efficient scaling:
    * - use theoretical scaling factors given by Mendel and Schwob ("Fast C-K-R Partitions of Sparse Graphs")
    * - use the precomputed MST as preprocessing
    * - embrace LCA DS for better performance and Quotient Graph Generating
 */
class EfficientCKR{
public:

    std::vector<double> scales;
    MendelScaling::UltrametricTree ultrametric;

    double diameter = 0.0;
    bool debug = false;

    std::vector<CKRLevel> m_levels;          // from finest (0) upward

    std::shared_ptr<TreeNode> getTree();

    EfficientCKR(IGraph& _g): g(_g) {
    }

    void updateEdgeDistances(const std::vector<double>& distances);
private:
    IGraph& g;
    std::vector<int> build_ckr_level(const IGraph& g, double Delta, CKRLevel& L);
    void preprocess();

    void computeScales();
    void computeLevelPartition(
        MendelScaling::QuotientLevel& Q,
        double Delta,
        CKRLevel& L);
    void buildTreeLevel(
        std::vector<std::shared_ptr<TreeNode>>& prev_nodes,
        const MendelScaling::QuotientLevel &Q,
        const CKRLevel& L,
        const double Delta);
    void finishTree(std::shared_ptr<TreeNode>& root, const std::vector<std::shared_ptr<TreeNode>>& prev_nodes);

};


#endif //OBLIVIOUSROUTING_EFFICIENT_CKR_H