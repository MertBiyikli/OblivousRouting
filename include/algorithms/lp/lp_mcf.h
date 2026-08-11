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
#include "core/errors.h"

using namespace operations_research;



class CMMF_Solver: public LP, public IOfflineSolver{
private:
    std::unordered_map<std::pair<int, int>, std::unordered_map<int,  MPVariable*>, PairHash> map_commodities2edge;
public:

    CMMF_Solver(optimized::Graph<EdgeData>& graph) : IOfflineSolver(graph), LP(graph.getNumNodes()) {
    }

    virtual Result<void> computeBasisFlows(AllPairRoutingTable& table) override;

    virtual Result<void> CreateVariables() override;

    virtual void CreateConstraints() override;
    virtual void SetObjective() override;
    virtual void storeFlow(AllPairRoutingTable& table) override;

    void PrintSolution();
    void AddDemandMap(const demands& d_map);

    Result<double> getCongestionForPassedDemandMap() const;
};


inline double computeRoutingSchemeCongestion(optimized::Graph<EdgeData>& _g,
                                             const std::unique_ptr<RoutingScheme>& routing_scheme,
                                             const demands& demand_map) {
    std::vector<double> congestion_per_edge(_g.getNumDirectedEdges(), 0.0);
    routing_scheme->routeDemands(congestion_per_edge, demand_map);
    double max_cong = routing_scheme->getMaxCongestion(congestion_per_edge);
    for (const auto& cong : congestion_per_edge)
        if (cong > max_cong) max_cong = cong;
    return max_cong;
}

inline Result<double> computeOfflineOptimalCongestion(optimized::Graph<EdgeData>& _g, const demands& demand_map) {
    CMMF_Solver mccf(_g);
    mccf.AddDemandMap(demand_map);
    auto offline_scheme = mccf.solve();
    auto congestion = mccf.getCongestionForPassedDemandMap();
    if (congestion) {
        return congestion.value();
    }else {
        return getError(congestion);
    }
}
#endif //OBLIVIOUSROUTING_LP_MCF_H