//
// Created by Mert on 27.10.26.
//

#ifndef E_ROUTING_COMMODITY_FLOW_DIAGNOSTICS_H
#define E_ROUTING_COMMODITY_FLOW_DIAGNOSTICS_H
#include "analysis/scenario_results.h"
#include "utils/demands.h"

struct CommodityFlowChange {
    int source;
    int target;
    double demand;
    double baseline_flow;
    double modified_flow;
    double baseline_load;
    double modified_load;
    double delta;
};

inline void analyzeCommodityFlows(const RoutingScheme& baseline, const RoutingScheme& modified, const demands& fixed_demands, int edge_id, std::size_t top_k = 10) {
    std::vector<CommodityFlowChange> changes;
    changes.reserve(fixed_demands.size());

    double baseline_total = 0.0;
    double modified_total = 0.0;

    for (std::size_t i = 0; i < fixed_demands.size(); ++i) {

        const auto [s, t] = fixed_demands.getDemandPair(i);

        const double demand = fixed_demands.getDemandValue(i);

        const double before = baseline.getFlow(edge_id, s, t);

        const double after = modified.getFlow(edge_id, s, t);

        // Diagnostic assumption:
        // signed coefficient magnitude represents
        // physical-link flow per demand unit.
        const double before_load = std::abs(before) * demand;

        const double after_load = std::abs(after) * demand;

        baseline_total += before_load;
        modified_total += after_load;

        changes.push_back({s, t, demand, before, after, before_load, after_load, after_load - before_load});
    }

    std::sort(changes.begin(), changes.end(), [](const auto& a, const auto& b) { return std::abs(a.delta) > std::abs(b.delta); });

    std::cout << "\n[COMMODITY FLOW DIAGNOSTICS]\n";
    std::cout << "Edge: " << edge_id << '\n';

    std::cout << "Baseline reconstructed load: " << baseline_total << '\n';

    std::cout << "Modified reconstructed load: " << modified_total << '\n';

    std::cout << "Total load difference: " << modified_total - baseline_total << '\n';

    std::cout << "\nTop commodity changes:\n";

    const auto count = std::min(top_k, changes.size());

    for (std::size_t i = 0; i < count; ++i) {
        const auto& c = changes[i];

        std::cout << c.source << " -> " << c.target << " | Demand: " << c.demand << " | Baseline coefficient: " << c.baseline_flow << " | Modified coefficient: " << c.modified_flow
                  << " | Baseline load: " << c.baseline_load << " | Modified load: " << c.modified_load << " | Delta: " << c.delta << '\n';
    }


}
#endif // E_ROUTING_COMMODITY_FLOW_DIAGNOSTICS_H
