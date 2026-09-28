#pragma once

#include "routing_result.h"
#include "core/config.h"
#include "core/errors.h"
#include "data_structures/graph/Igraph.h"
#include "data_structures/graph/graph.h"

class RoutingValidation {
public:
    [[nodiscard]] Result<void> input(const optimized::Graph<EdgeData>& _graph, const Config& _config ) const ;
    [[nodiscard]] Result<void> output(const IRoutingResult& _routing_results) const;

private:
    // input
    [[nodiscard]] Result<void> validateGraph(const optimized::Graph<EdgeData>& _graph) const;
    [[nodiscard]] Result<void> validateConfig(const Config& _config) const;
    [[nodiscard]] Result<void> validateDemand(const demands& _demands, const optimized::Graph<EdgeData>& graph) const;

    // output
    [[nodiscard]] Result<void> validateStatus(const IRoutingResult& _routing_results) const;
    [[nodiscard]] Result<void> validateTimeStats(const IRoutingResult& _routing_results) const;
    [[nodiscard]] Result<void> validateRouting(const IRoutingResult& _routing_results) const;
};

