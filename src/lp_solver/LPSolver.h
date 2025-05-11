//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_LPSOLVER_H
#define OBLIVOUSROUTING_LPSOLVER_H

#include "../solver/solver.h"

class LPSolver : public OblviviousRoutingSolver {
public:
    LPSolver() = default;

    void solve(const Graph& graph) override {
        // Implement the LP solver logic here
        // This is a placeholder implementation
        m_routingTable.resize(graph.numNodes(), std::vector<double>(graph.numNodes(), 0.0));
    }
};

#endif //OBLIVOUSROUTING_LPSOLVER_H
