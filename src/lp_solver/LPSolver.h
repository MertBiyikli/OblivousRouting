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
#include "../tree_based/frt/raecke_frt_transform.h"

/*
 * This is an attempt of implementing the LP solver for the oblivious routing problem.
 * Using the paper "Making Intra-Domain Routing Robust to Changing and
 * Uncertain Traffic Demands:
 * Understanding Fundamental Tradeoffs" from Applegate & Cohen
 */


class LPSolver : public LP {
private:
    std::unordered_map<std::tuple<int, int, int>, MPVariable*> p_e_ij;
    std::unordered_map< std::pair<int, int>, MPVariable*, PairHash> π_e_f;
    std::unordered_map< std::tuple< int , std::pair<int, int> > , MPVariable*> m_var_f_e_; // st;
    double max_cong = 0;

    virtual void CreateVariables(const GraphCSR& graph) override;
    virtual void CreateConstraints(const GraphCSR& graph) override;
    virtual void SetObjective() override;

    std::vector<int> sourceVertices;

public:
    //std::unordered_map<std::pair<int, int>, std::unordered_map<std::pair<int, int>, double>> f_st_e;
    LPSolver(){};


    double getMaximumCongestion(const GraphCSR& graph) const;
    bool SolveOptimalObliviousRouting(const GraphCSR& graph, bool AGGREGATE_CONGESTION = true);

    void PrintSolution(const GraphCSR& graph) override;
    void GetRoutingTable(const GraphCSR& graph);

    virtual void storeFlow() override;
    double getCongestion(DemandMap& demands, GraphCSR& g) const;

    void PrintCommoditiesPerEdge(const GraphCSR& graph);

    double getDemandWeightedCongestion(const GraphCSR& graph,
                                                 const DemandMap& demands);
};

#endif //OBLIVOUSROUTING_LPSOLVER_H
