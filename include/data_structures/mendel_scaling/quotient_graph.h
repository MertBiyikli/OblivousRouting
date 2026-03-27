//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_QUOTIENT_GRAPH_H
#define OBLIVIOUSROUTING_QUOTIENT_GRAPH_H
#include <cassert>
#include <memory>
#include <vector>

#include "ultrametric_tree.h"
#include "../graph/Igraph.h"

class QuotientLevel {
public:
    std::unique_ptr<IGraph> Gq;                                 // quotient graph at Δ
    std::vector<int> sigma_compact_of_v;      //original vertex v -> compact qid in [0..k-1]
    std::vector<std::vector<int>> members_of_q; // list of original vertices in each quotient node
};

/**
* For a given ultrametric tree and a scale Δ, we construct a quotient graph Gq where:
* - Each node in Gq corresponds to a cluster of original vertices that share the same ancestor at level Δ in the ultrametric tree.
* - There is an edge between two nodes in Gq if there is at least one edge in the original graph connecting any vertex in one cluster to any vertex in the other cluster, and the weight of that edge is within the sliding window [Δ/2n, Δ].
* - The capacity of the edge in Gq is the minimum capacity among all such edges in the original graph that connect the two clusters and fall within the sliding window.
*/
class QuotientConstruction {
public:
    int original_n = 0;
    std::vector<std::pair<int, double>> weights;

    void preprocessEdges(const IGraph& G);

    QuotientLevel constructQuotientGraph(
            const UltrametricTree& ultra,
            double Delta,
            IGraph& G);
};

#endif //OBLIVIOUSROUTING_QUOTIENT_GRAPH_H