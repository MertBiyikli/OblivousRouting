//
// Created by Mert on 22.10.26.
//

#ifndef E_ROUTING_FIXED_ROUTING_SCENARIO_H
#define E_ROUTING_FIXED_ROUTING_SCENARIO_H

#include <unordered_map>

#include "core/errors.h"
#include "visualization/visualization_result.h"

Result<RoutingVisualizationResult> transformScenario(const RoutingVisualizationResult& baseline, double demandScale, const std::unordered_map<int, double>& capacityFactors);

#endif // E_ROUTING_FIXED_ROUTING_SCENARIO_H
