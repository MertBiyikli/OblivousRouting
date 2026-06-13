//
// Created by Mert Biyikli on 12.06.26.
//

#ifndef OBLIVIOUSROUTING_OR_TOOLS_OPTIMIZER_H
#define OBLIVIOUSROUTING_OR_TOOLS_OPTIMIZER_H
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "ortools/linear_solver/linear_solver.h"
#include "load_optimizer.h"
#include "../../utils/demands.h"

class OrToolsSemiObliviousLoadOptimizer
    : public ISemiObliviousRoutingLoadOptimizer {
public:
    SemiObliviousOptimizationResult optimize(
        const IGraph& graph,
        const CandidateRoutingScheme& candidateScheme,
        const demands& demand
    ) override;
};

#endif //OBLIVIOUSROUTING_OR_TOOLS_OPTIMIZER_H