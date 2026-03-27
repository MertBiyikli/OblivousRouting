//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_LP_BASE_H
#define OBLIVIOUSROUTING_LP_BASE_H


#include "../solver.h"
#include "ortools/linear_solver/linear_solver.h"
#include <vector>
#include "../../data_structures/graph/Igraph.h"


using namespace operations_research;

class LP : public AllPairObliviousSolverBase {
public:
    bool debug = false;
    int n;
    MPSolver::ResultStatus status;
    std::unique_ptr<MPSolver> solver;
    MPVariable* alpha;
    std::vector<std::pair<int, int>> m_demands;

    LP(IGraph& graph):AllPairObliviousSolverBase(graph), solver(nullptr), alpha(nullptr) {
        init();
    }
    virtual ~LP() = default;
    void computeBasisFlows(AllPairRoutingTable &table) override;

    virtual void CreateVariables() = 0;
    virtual void CreateConstraints() = 0;
    virtual void SetObjective() = 0;
    virtual void PrintSolution() = 0;
    virtual void storeFlow(AllPairRoutingTable& table) = 0;


    void init();
    bool Run( AllPairRoutingTable& table);



};

#endif //OBLIVIOUSROUTING_LP_BASE_H