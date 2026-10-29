//
// Created by Mert on 28.10.26.
//

#ifndef E_ROUTING_SCENARIO_TRAFFIC_GROWTH_H
#define E_ROUTING_SCENARIO_TRAFFIC_GROWTH_H

#include "analysis/scenario_analysis.h"
#include "demand_vector.h"

class TrafficGrowth : public Scenario {
  public:
    TrafficGrowth() : Scenario() {}
    Result<void> verify(const optimized::Graph<EdgeData>& graph) const override;
    ScenarioKind kind() const noexcept override { return ScenarioKind::TRAFFIC_GROWTH; }

    Result<void> configure(const Config& cfg) override;

    Result<void> apply(optimized::Graph<EdgeData>& graph) const {
        return {}; // No network changes, only chaning traffic
    }

    Result<void> setGlobalScale(double scale);
    Result<void> growPerCommodity(int s, int t, double scale);
    Result<DemandVector> transform(const DemandVector& baseline) const;

private:

    double globalScale  = 1.0;
    std::unordered_map<std::pair<int, int>, double, PairHash> commodities;
};

#endif // E_ROUTING_SCENARIO_TRAFFIC_GROWTH_H
