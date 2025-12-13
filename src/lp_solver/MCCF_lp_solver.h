//
// Created by Mert Biyikli on 30.06.25.
//

#ifndef OBLIVIOUSROUTING_MCCF_LP_SOLVER_H
#define OBLIVIOUSROUTING_MCCF_LP_SOLVER_H


#include "../solver/solver.h"
#include "ortools/linear_solver/linear_solver.h"
#include <map>
#include <unordered_map>
#include <vector>
#include <tuple>
#include "../utils/hash.h"
#include "LP_Base.h"

/*
 * This class compute the Maximum Concurrent Flow Problem for a given set of demands.
 */

using namespace operations_research;

class CMMF_Solver: public LP{
private:
    std::unordered_map<std::pair<int, int>, double, PairHash> m_demands; // Flow variables for edges
    std::unordered_map<std::pair<int, int>, std::unordered_map<int,  MPVariable*>, PairHash> map_vertex2edge;
public:

    CMMF_Solver(IGraph& graph) : LP(graph) {
    }
    virtual void CreateVariables(const IGraph& graph) override;
    virtual void CreateConstraints(const IGraph& graph) override;
    virtual void SetObjective() override;
    void PrintSolution(const IGraph& graph) override;

    void AddDemands(const std::pair<int, int>& d, double value); // TODO: use here the Demand struct as defined in LPSolver.h
    virtual void storeFlow(EfficientRoutingTable& table) override;

};

#endif //OBLIVIOUSROUTING_MCCF_LP_SOLVER_H
