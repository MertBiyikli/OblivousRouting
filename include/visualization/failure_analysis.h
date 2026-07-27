//
// Created by Mert Biyikli on 24.07.26.
//

#ifndef OBLIVIOUSROUTING_EDGE_LOAD_ANALYZER_H
#define OBLIVIOUSROUTING_EDGE_LOAD_ANALYZER_H

#include <cstddef>
#include <string>
#include <vector>

#include "core/errors.h"
#include "data_structures/graph/Igraph.h"
#include "routing/routing_table.h"
#include "utils/demands.h"

struct AffectedCommodity {
    int source = -1;
    int target = -1;

    double demand = 0.0;

    /*
     * Fraction of the unit routing for this commodity that crosses the
     * failed undirected link.
     */
    double failed_fraction = 0.0;

    /*
     * demand × failed_fraction
     */
    double lost_traffic = 0.0;
};

struct LinkFailureSummary {
    int failed_edge_id = -1;
    int anti_edge_id = -1;

    int source = -1;
    int target = -1;

    double capacity = 0.0;

    std::size_t affected_commodities = 0;

    double affected_demand_volume = 0.0;
    double lost_traffic = 0.0;
    double total_demand = 0.0;

    double affected_demand_fraction = 0.0;
    double lost_traffic_fraction = 0.0;

    std::vector<AffectedCommodity> commodities;
};

class LinkFailureAnalyzer {
public:

    static Result<LinkFailureSummary> analyze(const IGraph& graph,const RoutingScheme& scheme,const demands& demand_map,int failed_edge);
};

#endif //OBLIVIOUSROUTING_EDGE_LOAD_ANALYZER_H