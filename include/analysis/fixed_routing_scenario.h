//
// Created by Mert on 22.10.26.
//

#ifndef E_ROUTING_FIXED_ROUTING_SCENARIO_H
#define E_ROUTING_FIXED_ROUTING_SCENARIO_H

#include "core/errors.h"
#include "scenario_results.h"
#include <unordered_map>

Result<NetworkData> transformScenario(const NetworkData& baseline, double demandScale, const std::unordered_map<int, double>& capacityFactors);

#endif // E_ROUTING_FIXED_ROUTING_SCENARIO_H
