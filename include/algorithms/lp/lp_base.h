//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_LP_BASE_H
#define OBLIVIOUSROUTING_LP_BASE_H

#include "ortools/linear_solver/linear_solver.h"
#include "routing/storage/allpair_routing_table.h"
#include "core/errors.h"
#include <vector>

using namespace operations_research;

class LP {
public:
    bool debug = false;
    int n;
    MPSolver::ResultStatus status;
    std::unique_ptr<MPSolver> solver;
    MPVariable* alpha;
    std::vector<std::pair<int, int>> m_demands;

    LP(const int& _n): solver(nullptr), alpha(nullptr)  {
        n = _n;
    }

    virtual ~LP() = default;

    Result<void> initSolver() {
        solver = std::unique_ptr<MPSolver>(MPSolver::CreateSolver("GLOP"));
        if (!solver) {
            return makeErrorMessage(ErrorCode::InvalidSolver, "LP Solver failed to initializes.");
        }else {
            return {};
        }
    }
    virtual Result<void> CreateVariables() = 0;
    virtual void CreateConstraints() = 0;
    virtual void SetObjective() = 0;
    virtual void storeFlow(AllPairRoutingTable& table) = 0;

    bool Run( AllPairRoutingTable& table);



};

#endif //OBLIVIOUSROUTING_LP_BASE_H