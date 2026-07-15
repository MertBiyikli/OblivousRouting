#pragma once

#include "../hierarchy_results.h"

#include "algorithms/oblivious/mwu/oracle/electrical/laplacian_solver.h"
#include "core/errors.h"
#include "data_structures/graph/Igraph.h"

#include <vector>

struct LocalElectricalFlowResult {
    // Indexed like cluster.induced_edges.
    std::vector<double> edge_flow;

    double max_congestion = 0.0;
    double max_conservation_error = 0.0;
};

class LocalElectricalSolver {
public:
    /*
     * local_imbalance[i] corresponds to
     * cluster.original_vertices[i].
     */

    Result<LocalElectricalFlowResult> routeDemand(const IGraph& graph,const HierarchyCluster& cluster,const std::vector<double>& local_imbalance) const;
};