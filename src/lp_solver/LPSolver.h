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
    std::unordered_map< std::pair<int, int>, MPVariable*, PairHash> π_e_f;
    std::unordered_map< std::tuple< int , std::pair<int, int> > , MPVariable*> m_var_f_e_; // st;
    double max_cong = 0;

    virtual void CreateVariables(const IGraph& graph) override;
    virtual void CreateConstraints(const IGraph& graph) override;
    virtual void SetObjective() override;

    std::vector<int> sourceVertices;

public:
    //std::unordered_map<std::pair<int, int>, std::unordered_map<std::pair<int, int>, double>> f_st_e;

    LPSolver(IGraph& graph):LP(graph){};


    double getMaximumCongestion(const IGraph& graph) const;
    bool SolveOptimalObliviousRouting(const IGraph& graph, bool AGGREGATE_CONGESTION = true);

    void PrintSolution(const IGraph& graph) override;
    void GetRoutingTable(const IGraph& graph);

    virtual void storeFlow(AllPairRoutingTable& table) override;
    double getCongestion(DemandMap& demands, IGraph& g) const;

    void PrintCommoditiesPerEdge(const IGraph& graph);

    double getDemandWeightedCongestion(const IGraph& graph,
                                                 const DemandMap& demands);
};

#endif //OBLIVOUSROUTING_LPSOLVER_H
