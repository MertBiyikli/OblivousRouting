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
#include "visualization/failure_analysis.h"

struct VisualizationNode {
    int id = -1;
    std::string label;
};

struct VisualizationEdge {
    int id = -1;
    int source = -1;
    int target = -1;

    double capacity = 0.0;
    double load = 0.0;
    double utilization = 0.0;

    bool enabled = true;
};


struct RoutingAnalysisSummary {
    double maximum_congestion = 0.0;
    double average_utilization = 0.0;
    double total_routed_demand = 0.0;

    std::size_t overloaded_edges = 0;
};

struct RoutingVisualizationResult {
    std::string graph_name;
    std::string solver_name;
    std::string demand_model;

    std::vector<VisualizationNode> nodes;
    std::vector<VisualizationEdge> edges;

    RoutingAnalysisSummary summary;
    /*
     * One precomputed static failure analysis per physical edge.
     */
    std::vector<LinkFailureSummary> link_failures;
};

class RoutingAnalyzer {
public:
    static Result<RoutingVisualizationResult> analyze(const optimized::Graph<EdgeData>& graph,const RoutingScheme& scheme,const demands& demand_map,std::string graph_name,std::string solver_name,std::string demand_model);
    static Result<void> analyzeSingleLinkFailures(const optimized::Graph<EdgeData>& graph,const RoutingScheme& scheme,const demands& demand_map,RoutingVisualizationResult& result);
};
#endif //OBLIVIOUSROUTING_VISUALIZATION_RESULT_H