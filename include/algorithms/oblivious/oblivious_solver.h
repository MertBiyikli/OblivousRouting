//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_SOLVER_H
#define OBLIVOUSROUTING_SOLVER_H

#include "routing/routing_table.h"
#include "routing/storage/linear_routing_table.h"
#include "routing/storage/allpair_routing_table.h"
#include "core/solver.h"
#include "data_structures/graph/Igraph.h"
#include <memory>


class ILinearObliviousSolverBase : public ISolver {
protected:
    int root = 0;

public:
    ILinearObliviousSolverBase(IGraph& g, int root)
        : ISolver(g), root(root) {}

    std::unique_ptr<RoutingScheme> solve() final {
        LinearRoutingTable table;
        table.init(graph);

        computeBasisFlows(table);

        graph.resetEdgeDistance();

        return std::make_unique<LinearRoutingScheme>(
            graph,
            root,
            std::move(table)
        );
    }

    int getRootNode() const {
        return root;
    }

protected:
    virtual void computeBasisFlows(LinearRoutingTable& table) = 0;
};


class IAllPairObliviousSolverBase : public ISolver {
public:
    explicit IAllPairObliviousSolverBase(IGraph& g)
        : ISolver(g) {}

    std::unique_ptr<RoutingScheme> solve() {
        AllPairRoutingTable table;
        table.init(graph);

        computeBasisFlows(table);

        graph.resetEdgeDistance();

        return std::make_unique<AllPairRoutingScheme>(
            graph,
            std::move(table)
        );
    }

protected:
    virtual void computeBasisFlows(AllPairRoutingTable& table) = 0;
};



#endif //OBLIVOUSROUTING_SOLVER_H
