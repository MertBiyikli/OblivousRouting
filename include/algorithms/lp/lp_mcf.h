//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_LP_MCF_H
#define OBLIVIOUSROUTING_LP_MCF_H

#include "lp_base.h"
#include "ortools/linear_solver/linear_solver.h"
#include "../../utils/hash.h"
#include "../../utils/demands.h"

using namespace operations_research;

class CMMF_Solver: public LP{
private:
    std::unordered_map<std::pair<int, int>, double, PairHash> m_demand_map; // Flow variables for edges
    std::unordered_map<std::pair<int, int>, std::unordered_map<int,  MPVariable*>, PairHash> map_vertex2edge;
public:

    CMMF_Solver(IGraph& graph) : LP(graph) {
    }
    virtual void CreateVariables() override;
    virtual void CreateConstraints() override;
    virtual void SetObjective() override;
    void PrintSolution() override;

    void AddDemandMap(const demands& d_map);
    void AddDemands(const std::pair<int, int>& d, double value);
    virtual void storeFlow(AllPairRoutingTable& table) override;

    double getCongestionForPassedDemandMap();
};

#endif //OBLIVIOUSROUTING_LP_MCF_H