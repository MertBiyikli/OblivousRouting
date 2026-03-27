//
// Created by Mert Biyikli on 17.09.25.
//

#ifndef OBLIVIOUSROUTING_MST_ALGO_H
#define OBLIVIOUSROUTING_MST_ALGO_H

#include <vector>
#include <tuple>
#include "../../../../../data_structures/graph/Igraph.h"
#include "../../../../../data_structures/hst/pointer_hst.h"
#include "../../../../../data_structures/hst/flat_hst.h"
#include "../../../../../data_structures/union_find/union_find.h"


/**
 * compute minimum spanning tree using Kurskal's algorithm
 */
class MST {
    int n;
    std::vector<std::tuple<double,int,int>> keyed;


public:
    MST() = delete;
    explicit MST(IGraph& g);
    // Build a random MST edge set using Kruskal with random priorities
    std::vector<std::pair<int,int>> computeMST();

};

#endif // OBLIVIOUSROUTING_MST_ALGO_H
