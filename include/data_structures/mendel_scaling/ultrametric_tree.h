//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_ULTRAMETRIC_TREE_H
#define OBLIVIOUSROUTING_ULTRAMETRIC_TREE_H
#include <vector>

// Each original vertex v is a leaf at index v (0..n-1).
// Internal nodes are appended at indices n .. n+merges-1.
struct UltraMetricNode {
    int parent = -1;
    std::vector<int> child;
    double Gamma = 0.0;   // label Γ(u) (nondecreasing toward root)
};


/**
 * Ultrametric tree structure for Mendel-Naor scaling.
 * It is built from an MST of the original graph, where each union operation in the MST construction
 * corresponds to creating a new internal node in the ultrametric tree.
 * The Γ label of each internal node is set to the weight of the edge that caused the union.
 *
 * Importantly, note that in the original paper, they introduce a factor (n-1) to multiply the weights with
 * to ensure the ultrametric property holds. However, in our implementation, we can directly use the original weights
 * since in the original paper they first approximate the graph by a spanner which introduce the additional (n-1) approximation factor.
 */
class UltrametricTree {
    public:
    int n = 0;                // number of original nodes
    int N = 0;                // total nodes = leaves (original node) + internal
    int root = -1;
    std::vector<UltraMetricNode> T;

    // Binary lifting tables for level-ancestor-like jumps
    std::vector<std::vector<int>> up;     // up[k][u] = 2^k-th ancestor of u (or -1)
    std::vector<double> gamma;            // Γ(u), cached for easy access
    std::vector<int> depth;               // depth in the ultrametric tree


    /** Build ultrametric from an MST using Kruskal-like dendrogram construction:
    * Sort MST edges by weight ascending; each union creates a new internal node
    * with Γ = edge weight*; connect components as its children.
    */
    void buildFromMST(int n_, const std::vector<std::pair<int, int>>& mst_edges, const std::vector<double>& mst_w);

    void preprocessLifting();

    /** Return the highest ancestor of leaf v whose Γ <= threshold.
    * If leaf’s Γ(leaf)=0 is already > threshold (never happens), it returns the leaf itself.
    */
    int sigmaDelta(int v_leaf, double Delta) const;
};

#endif //OBLIVIOUSROUTING_ULTRAMETRIC_TREE_H