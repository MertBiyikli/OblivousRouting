//
// Created by Mert Biyikli on 15.08.25.
//

#ifndef OBLIVIOUSROUTING_LP_OBLIVIOUS_RATIO_H
#define OBLIVIOUSROUTING_LP_OBLIVIOUS_RATIO_H

#include "ortools/linear_solver/linear_solver.h"
#include <map>
#include <vector>
#include <tuple>
#include "../../src/utils/hash.h"
#include "../../src/graph.h"
#include "../../src/lp_solver/LP_Base.h"



using namespace operations_research;

class ObliviousRatio{
    bool debug = false;
    std::unordered_set< std::pair<int, int> > demands;
    std::vector< std::pair<int, int> > edges;
    Graph g;
    // This is the oblivious routing strategy for which the oblivious ratio should be computed for
    std::unordered_map<std::pair<int,int>, std::unordered_map<std::pair<int,int>, double>> routing;

    // Worst case LP result
    double max_congestion = 0;
    std::unordered_map<std::pair<int, int>, double> demand_vector;

    public:
    void init(Graph& g, const std::unordered_map<std::pair<int,int>, std::unordered_map<std::pair<int,int>, double>>& routing);
    double solve();
    double solveWorstCaseDemandLPPerEdge(const std::pair<int, int>& edge);
};

#endif //OBLIVIOUSROUTING_LP_OBLIVIOUS_RATIO_H