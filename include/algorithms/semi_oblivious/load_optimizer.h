//
// Created by Mert Biyikli on 12.06.26.
//

#ifndef OBLIVIOUSROUTING_LOAD_OPTIMIZER_H
#define OBLIVIOUSROUTING_LOAD_OPTIMIZER_H

#include "../../data_structures/graph/Igraph.h"
#include "../../utils/demands.h"
#include "../../core/routing_scheme.h"
#include "candidate_routing_scheme.h"
#include <optional>


struct SemiObliviousOptimizationResult {
    std::unique_ptr<RoutingScheme> scheme;
    double lambda = -1.0;
};

class ISemiObliviousRoutingLoadOptimizer {
public:
    virtual ~ISemiObliviousRoutingLoadOptimizer() = default;
    virtual SemiObliviousOptimizationResult optimize(
        const IGraph& g,
        const CandidateRoutingScheme& candidate_routing_scheme,
        const demands& Demands) = 0;
};

#endif //OBLIVIOUSROUTING_LOAD_OPTIMIZER_H