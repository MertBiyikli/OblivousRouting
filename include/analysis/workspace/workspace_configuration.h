//
// Created by Mert on 30.10.26.
//

#ifndef E_ROUTING_WORKSPACE_CONFIGURATION_H
#define E_ROUTING_WORKSPACE_CONFIGURATION_H

#include "core/config.h"
#include <vector>

struct CommodityModification {
    int source = -1;
    int target = -1;
    double factor = 1.0;
};

struct WorkspaceConfiguration {
    std::vector<_CapacityUpgrade> capacity_upgrades;

    double demand_scale = 1.0;

    std::vector<CommodityModification> commodity_overrides;

    // Reset all user-defined modifications.
    void reset() {
        capacity_upgrades.clear();
        demand_scale = 1.0;
        commodity_overrides.clear();
    }

    bool isModified() const noexcept { return !capacity_upgrades.empty() || demand_scale != 1.0 || !commodity_overrides.empty(); }
};


#endif // E_ROUTING_WORKSPACE_CONFIGURATION_H
