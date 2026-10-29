//
// Created by Mert on 28.10.26.
//

#ifndef E_ROUTING_SCENARIO_CAPACITY_UPGRADE_H
#define E_ROUTING_SCENARIO_CAPACITY_UPGRADE_H
#include "analysis/scenario_analysis.h"



class CapacityUpgrade : public Scenario {
public:
    CapacityUpgrade() : Scenario() {}


    Result<void> verify(const optimized::Graph<EdgeData>& graph) const override;
    Result<void> apply(optimized::Graph<EdgeData>& graph) const override;
    ScenarioKind kind() const noexcept override { return ScenarioKind::CAPACITY_UPGRADE; }
    Result<void> configure(const Config& cfg) override;

    Result<void> upgrade(const Config& cfg);
    Result<void> upgrade(int edge_id, double scale);

    std::unordered_map<int, double> upgrades;
};

#endif // E_ROUTING_SCENARIO_CAPACITY_UPGRADE_H
