//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_SOLVER_H
#define OBLIVOUSROUTING_SOLVER_H

#include <vector>
#include "../graph.h"
#include "../utils/hash.h"

class ObliviousRoutingSolver { // TODO: rename the class name
public:
    std::unordered_map<std::pair<int, int> , std::unordered_map<std::pair<int, int>, double >> f_e_st;
    std::vector<std::vector<double>> m_routingTable;
    ObliviousRoutingSolver() = default;
    virtual ~ObliviousRoutingSolver() = default;

    // ToDo: Think of an efficient way of storing the routing table
    virtual void solve(const Graph& graph) = 0;
    virtual void storeFlow() = 0;
};

#endif //OBLIVOUSROUTING_SOLVER_H
