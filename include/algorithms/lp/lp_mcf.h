//
// Created by Mert Biyikli on 25.03.26.
//

#ifndef OBLIVIOUSROUTING_LP_MCF_H
#define OBLIVIOUSROUTING_LP_MCF_H

#include "lp_base.h"
#include "offline_solver.h"
#include "ortools/linear_solver/linear_solver.h"
#include "../../utils/hash.h"
#include "../../utils/demands.h"

using namespace operations_research;

class CMMF_Solver: public LP, public IOfflineSolver{
private:
    std::unordered_map<std::pair<int, int>, std::unordered_map<int,  MPVariable*>, PairHash> map_commodities2edge;
public:

    CMMF_Solver(IGraph& graph) : IOfflineSolver(graph) {}

    virtual void computeBasisFlows(AllPairRoutingTable& table) override;

    virtual void CreateVariables() override;
    virtual void CreateConstraints() override;
    virtual void SetObjective() override;
    virtual void storeFlow(AllPairRoutingTable& table) override;

    void PrintSolution();
    void AddDemandMap(const demands& d_map);
    void AddDemands(const std::pair<int, int>& d, double value);

    double getCongestionForPassedDemandMap() const;
};

#endif //OBLIVIOUSROUTING_LP_MCF_H