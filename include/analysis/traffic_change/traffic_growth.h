//
// Created by Mert on 28.10.26.
//

#ifndef E_ROUTING_SCENARIO_TRAFFIC_GROWTH_H
#define E_ROUTING_SCENARIO_TRAFFIC_GROWTH_H

#include "analysis/scenario_analysis.h"
/*
class TrafficGrowth : public Scenario {
public:
    TrafficGrowth() : Scenario() {}
    Result<void> verify(const optimized::Graph<EdgeData>& graph) const override;
    ScenarioKind kind() const noexcept override { return ScenarioKind::TRAFFIC_GROWTH; }

    Result<void> setGlobalScale(double scale);
    Result<void> growPerCommodity(int s, int t, double scale);
    double globalScale{};
    std::unordered_map<std::pair<int, int>, double, PairHash> commodities;
};*/

#endif // E_ROUTING_SCENARIO_TRAFFIC_GROWTH_H
