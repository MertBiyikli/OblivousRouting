//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_LP_AC_H
#define OBLIVIOUSROUTING_LP_AC_H

#include "lp_base.h"
#include "ortools/linear_solver/linear_solver.h"
#include "../../utils/hash.h"
#include <vector>
#include <tuple>


/*
 * This is an implemention of the LP solver for the oblivious routing problem.
 * Using the paper "Making Intra-Domain Routing Robust to Changing and
 * Uncertain Traffic Demands:
 * Understanding Fundamental Tradeoffs" from Applegate & Cohen
 */
class LPSolver : public LP {
private:
    std::unordered_map<std::tuple<int, int, int>, MPVariable*, TrippleTuple> p_e_ij;
    std::unordered_map< std::pair<int, int>, MPVariable*, PairHash> π_e_f;
    std::unordered_map< std::tuple< int , std::pair<int, int> > , MPVariable*, TuplePair> m_var_f_e_; // st;
    double max_cong = 0;

    virtual void CreateVariables() override;
    virtual void CreateConstraints() override;
    virtual void SetObjective() override;

    std::vector<int> sourceVertices;

public:

    LPSolver(IGraph& graph):LP(graph){};

    virtual void storeFlow(AllPairRoutingTable& table) override;

    void PrintSolution() override;


};

#endif //OBLIVIOUSROUTING_LP_AC_H