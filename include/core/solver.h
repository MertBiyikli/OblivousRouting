//
// Created by Mert Biyikli on 22.06.26.
//

#ifndef OBLIVIOUSROUTING_SOLVER_H
#define OBLIVIOUSROUTING_SOLVER_H

#include "routing/routing_table.h"
#include "data_structures/graph/Igraph.h"
#include <memory>

class ISolver {
protected:
    IGraph& graph;

public:
    explicit ISolver(IGraph& g) : graph(g) {}
    virtual ~ISolver() = default;

    virtual std::unique_ptr<RoutingScheme> solve() = 0;
};
#endif //OBLIVIOUSROUTING_SOLVER_H