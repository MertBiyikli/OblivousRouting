//
// Created by Mert Biyikli on 25.03.26.
//
#include "../../../include/algorithms/lp/lp_base.h"

void LP::init() {
    n = this->graph.getNumNodes();

    if (!solver) {
        solver.reset(MPSolver::CreateSolver("GLOP"));
    }
    if (!solver) {
        throw std::runtime_error("GLOP solver unavailable.");
    }
}

void LP::computeBasisFlows(AllPairRoutingTable &table) {
    if (!Run( table)) {
        throw std::runtime_error("LP failed to find a valid solution.");
    }
}



bool LP::Run( AllPairRoutingTable& table) {
    CreateVariables();
    CreateConstraints();
    SetObjective();
    // === Solve the LP ===
    status = solver->Solve();
    if (status != MPSolver::OPTIMAL) {
        std::cerr << "WARNING: LP is infeasible or unbounded.\n";
        return false;
    } else {
        storeFlow(table);
        return true;
    }
}

