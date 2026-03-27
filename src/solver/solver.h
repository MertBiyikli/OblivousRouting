//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_SOLVER_H
#define OBLIVOUSROUTING_SOLVER_H

#include "routing_table.h"
#include <memory>


class ObliviousRoutingSolver {
public:
    virtual ~ObliviousRoutingSolver() = default;
    virtual std::unique_ptr<RoutingScheme> solve() = 0;
    virtual void setNumThreads(int num_threads) {};
};

class LinearObliviousSolverBase : public ObliviousRoutingSolver {
public:

    LinearObliviousSolverBase(IGraph& _g, int _root)
        : graph(_g), root(_root) {}

    std::unique_ptr<RoutingScheme> solve() override {
        LinearRoutingTable table;
        table.init(graph);

        // solver must fill table via computeBasisFlows
        computeBasisFlows(table);
        graph.resetEdgeDistance();

        return std::make_unique<LinearRoutingScheme>(
            graph,
            root,
            std::move(table));
    }

protected:
    IGraph& graph;
    int root;

    virtual void computeBasisFlows(LinearRoutingTable& table) = 0;
};

class AllPairObliviousSolverBase : public ObliviousRoutingSolver {
public:
    AllPairObliviousSolverBase(IGraph& _g):
        graph(_g) {}


    std::unique_ptr<RoutingScheme> solve() override {
        AllPairRoutingTable table;
        table.init(graph);

        computeBasisFlows(table);

        graph.resetEdgeDistance();

        return std::make_unique<AllPairRoutingScheme>(
            graph,
            std::move(table));
    }

protected:
    IGraph& graph;
    virtual void computeBasisFlows(AllPairRoutingTable& table) = 0;
};

#endif //OBLIVOUSROUTING_SOLVER_H
