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
 * This is an attempt of implementation of an LP solver for the oblivious routing problem.
 * Using the paper "Making Intra-Domain Routing Robust to Changing and
 * Uncertain Traffic Demands:
 * Understanding Fundamental Tradeoffs" from Applegate & Cohen
 */


class LPSolver : public OblviviousRoutingSolver, public LP {
private:
    std::unordered_map<std::tuple<int, int, int>, MPVariable*> p_e_ij;
    std::unordered_map< std::pair<int, int>, MPVariable*> π_e_f;
    std::unordered_map< std::tuple< int , std::pair<int, int> > , MPVariable*> f_e_st;

    double max_cong = 0;

    virtual void CreateVariables(const DiGraph& graph) override;
    virtual void CreateConstraints(const DiGraph& graph) override;

    std::vector<int> sourceVertices;

public:
    LPSolver():OblviviousRoutingSolver(){};
    void solve(const Graph& graph) override;


    bool SolveOptimalObliviousRouting(const DiGraph& graph, bool AGGREGATE_CONGESTION = true);

    void PrintObliviousRoutingSolution(const DiGraph& graph);
    void GetRoutingTable(const DiGraph& graph);
};

#endif //OBLIVOUSROUTING_LPSOLVER_H
