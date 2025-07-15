//
// Created by Mert Biyikli on 01.07.25.
//

#ifndef OBLIVIOUSROUTING_LP_BASE_H
#define OBLIVIOUSROUTING_LP_BASE_H

#include "../solver/solver.h"
#include "ortools/linear_solver/linear_solver.h"
#include <map>
#include <vector>
#include <tuple>
#include "../utils/hash.h"



using namespace operations_research;

class LP{
public:
    std::unique_ptr<MPSolver> solver;
    MPVariable* alpha;
    std::vector<Demand> demands;

    LP() : solver(nullptr), alpha(nullptr) {}
    virtual ~LP() = default;

    virtual void CreateVariables(const DiGraph& graph) = 0;
    virtual void CreateConstraints(const DiGraph& graph) = 0;

    void Run(const Graph& graph) {

        if (!solver) {
            solver.reset(MPSolver::CreateSolver("GLOP"));
        }
        if (!solver) {
            throw std::runtime_error("GLOP solver unavailable.");
        }
        auto& digraph = graph.GetDiGraph();
        CreateVariables(digraph);
        CreateConstraints(digraph);

        // === Objective: maximize alpha ===
        solver->MutableObjective()->SetCoefficient(alpha, 1);
        solver->MutableObjective()->SetMinimization();

        // === Solve the LP ===
        if (solver->Solve() != MPSolver::OPTIMAL) {
            std::cerr << "WARNING: LP is infeasible or unbounded.\n";
            return;
        } else {
            std::cout << "Optimal congestion alpha: " << alpha->solution_value() << std::endl;
        }

    }

};

#endif //OBLIVIOUSROUTING_LP_BASE_H
