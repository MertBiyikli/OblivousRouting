//
// Created by Mert Biyikli on 24.07.26.
//

#ifndef OBLIVIOUSROUTING_VISUALIZATION_RESULT_H
#define OBLIVIOUSROUTING_VISUALIZATION_RESULT_H


#include <cstddef>
#include <string>
#include <vector>

#include "data_structures/graph/Igraph.h"
#include "routing/routing_table.h"
#include "utils/demands.h"
#include "core/errors.h"
#include "routing/routing_result.h"
#include "visualization/failure_analysis.h"




class RoutingAnalyzer {
public:
    static Result<RoutingVisualizationResult> analyze(const optimized::Graph<EdgeData>& graph,const RoutingScheme& scheme,const demands& demand_map, const IRoutingResult& result,std::string demand_model);
    static Result<void> analyzeSingleLinkFailures(const optimized::Graph<EdgeData>& graph,const RoutingScheme& scheme,const demands& demand_map,RoutingVisualizationResult& result);
};
#endif //OBLIVIOUSROUTING_VISUALIZATION_RESULT_H