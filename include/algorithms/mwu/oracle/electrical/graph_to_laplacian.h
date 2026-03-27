//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_GRAPH_TO_LAPLACIAN_H
#define OBLIVIOUSROUTING_GRAPH_TO_LAPLACIAN_H

#include <vector>
#include "../../../../data_structures/graph/Igraph.h"

class GraphToLaplacian {
public:
    std::vector<int> to, from, head; // to[e], from[e] gives the endpoints of edge e
    std::vector<double> weights; // weights[e] gives the weight of edge e
    std::vector<int> laplacian_indices; // laplacian_indices[e] gives the index in the Laplacian matrix values array
    std::vector<int> laplacian_indices_for_diagonal_elements; // laplacian_indices_red[e] gives the index in the reduced Laplacian matrix values array

    void init(const IGraph& graph);

    double getEdgeWeight(int e) const;
    double getEdgeWeight(int u, int v) const;
    void setEdgeWeight(int e, double w);
    void setEdgeWeight(int u, int v, double w);
    void setLaplacianIndex(int u, int v, int idx);
    int getLaplacianIndex(int u, int v) const;
};
#endif //OBLIVIOUSROUTING_GRAPH_TO_LAPLACIAN_H