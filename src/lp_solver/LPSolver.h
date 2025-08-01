//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_LPSOLVER_H
#define OBLIVOUSROUTING_LPSOLVER_H

#include "../solver/solver.h"
#include "ortools/linear_solver/linear_solver.h"
#include <map>
#include <vector>
#include <tuple>
#include "../utils/hash.h"
#include "LP_Base.h"
/*
 * This is an attempt of implementing the LP solver for the oblivious routing problem.
 * Using the paper "Making Intra-Domain Routing Robust to Changing and
 * Uncertain Traffic Demands:
 * Understanding Fundamental Tradeoffs" from Applegate & Cohen
 */


class LPSolver : public LP {
private:
    std::unordered_map<std::tuple<int, int, int>, MPVariable*> p_e_ij;
    std::unordered_map< std::pair<int, int>, MPVariable*> π_e_f;
    std::unordered_map< std::tuple< int , std::pair<int, int> > , MPVariable*> f_e_st;

    double max_cong = 0;

    virtual void CreateVariables(const RaeckeGraph& graph) override;
    virtual void CreateConstraints(const RaeckeGraph& graph) override;
    virtual void SetObjective() override;

    std::vector<int> sourceVertices;

public:
    LPSolver(){};


    double getMaximumCongestion(const RaeckeGraph& graph) const;
    bool SolveOptimalObliviousRouting(const RaeckeGraph& graph, bool AGGREGATE_CONGESTION = true);

    void PrintSolution(const RaeckeGraph& graph) override;
    void GetRoutingTable(const RaeckeGraph& graph);

    void PrintCommoditiesPerEdge(const RaeckeGraph& graph);
};

#endif //OBLIVOUSROUTING_LPSOLVER_H
