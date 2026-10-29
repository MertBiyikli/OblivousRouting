//
// Created by Mert on 28.10.26.
//

#ifndef E_ROUTING_SCENARIO_LINK_FAILURE_H
#define E_ROUTING_SCENARIO_LINK_FAILURE_H
#include "analysis/scenario_analysis.h"

class LinkFailure : public Scenario {
public:
    LinkFailure() : Scenario() {}

    Result<void> verify(const optimized::Graph<EdgeData>& graph) const override;
    Result<void> apply(optimized::Graph<EdgeData>& g) const override;
    ScenarioKind kind() const noexcept override { return ScenarioKind::LINK_FAILURE; }
    Result<void> configure(const Config& cfg) override;
    Result<void> linkFailure(int edge_id);

    std::unordered_set<int> failed_links;
};



#endif // E_ROUTING_SCENARIO_LINK_FAILURE_H
