//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_SOLVER_H
#define OBLIVOUSROUTING_SOLVER_H

#include <vector>
#include "../graph.h"

class OblviviousRoutingSolver { // TODO: rename the class name
public:
    std::vector<std::vector<double>> m_routingTable;
    OblviviousRoutingSolver() = default;
    virtual ~OblviviousRoutingSolver() = default;

    // ToDo: Think of an efficient way of storing the routing table
    virtual void solve(const Graph& graph) = 0;
};

#endif //OBLIVOUSROUTING_SOLVER_H
