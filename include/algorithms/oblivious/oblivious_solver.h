//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_SOLVER_H
#define OBLIVOUSROUTING_SOLVER_H

#include "routing/routing_table.h"
#include "routing/storage/linear_routing_table.h"
#include "routing/storage/allpair_routing_table.h"
#include "core/solver.h"
#include "core/errors.h"
#include "data_structures/graph/Igraph.h"
#include <memory>


class ILinearObliviousSolverBase : public ISolver {
protected:
    int root = 0;

public:
    ILinearObliviousSolverBase(optimized::Graph<EdgeData>& g, int root)
        : ISolver(g), root(root) {}

    Result<std::unique_ptr<RoutingScheme>> solve() final {
        LinearRoutingTable table;
        table.init(graph);

        auto basis_flow = computeBasisFlows(table);
        if ( !basis_flow) {
            return getError(basis_flow);
        }

        resetEdgeDistance();

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
    virtual Result<void> computeBasisFlows(LinearRoutingTable& table) = 0;
};


class IAllPairObliviousSolverBase : public ISolver {
public:
    explicit IAllPairObliviousSolverBase(optimized::Graph<EdgeData>& g)
        : ISolver(g) {}

    Result<std::unique_ptr<RoutingScheme>> solve() {
        AllPairRoutingTable table;
        table.init(graph);

        auto basis_flow = computeBasisFlows(table);
        if ( !basis_flow) {
            return getError(basis_flow);
        }


        return std::make_unique<AllPairRoutingScheme>(
            graph,
            std::move(table)
        );
    }

protected:
    virtual Result<void> computeBasisFlows(AllPairRoutingTable& table) = 0;
};



#endif //OBLIVOUSROUTING_SOLVER_H
