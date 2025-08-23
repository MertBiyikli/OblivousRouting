//
// Created by Mert Biyikli on 18.08.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_SOLVER_H
#define OBLIVIOUSROUTING_RAECKE_SOLVER_H

#include "../solver/solver.h"
#include "raecke_tree_decomp.h"
#include "raecke_transform.h"

class RaeckeSolver : public ObliviousRoutingSolver {
    bool debug = false;
    RaeckeFRT raeckeFRT;
    RaeckeSolutionTransform raeckeTransform;
public:

    void solve(const Graph& graph) override;
    void storeFlow() override;
    double getCongestion();
};

#endif //OBLIVIOUSROUTING_RAECKE_SOLVER_H