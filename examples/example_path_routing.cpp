/**
 * Example: Using PathRoutingTable to track routing paths and evaluate congestion
 *
 * This example shows how to:
 * 1. Create a PathRoutingTable
 * 2. Add routing paths with flow
 * 3. Evaluate congestion metrics
 * 4. Reconstruct the linear flow representation
 */

#include "include/core/path_routing_table.h"
#include "include/core/routing_table.h"
#include "include/data_structures/graph/Igraph.h"

// Example usage (pseudocode):
/*

void example_path_routing() {
    // Load graph
    IGraph graph;
    graph.loadFromLGF("input.lgf");

    // Create path routing table
    PathRoutingTable path_table;
    path_table.init(graph);

    // Example: Add paths for source 1
    int source = 1;
    double flow_per_path = 0.5;

    // Path 1: 1 -> 2 -> 3 -> 0 (root)
    std::vector<int> path1 = {
        graph.getEdgeId(1, 2),
        graph.getEdgeId(2, 3),
        graph.getEdgeId(3, 0)
    };
    path_table.addPath(source, path1, flow_per_path);

    // Path 2: 1 -> 4 -> 0 (root)
    std::vector<int> path2 = {
        graph.getEdgeId(1, 4),
        graph.getEdgeId(4, 0)
    };
    path_table.addPath(source, path2, flow_per_path);

    // Similarly add paths for other sources...

    // Evaluate congestion
    double max_congestion = path_table.evaluateCongestion(graph);
    std::cout << "Max edge congestion: " << max_congestion << "\n";

    // Get detailed statistics
    auto stats = path_table.getCongestionStats(graph);
    std::cout << "Avg congestion: " << stats["avg"] << "\n";
    std::cout << "Std dev: " << stats["std_dev"] << "\n";

    // Get edge loads for all edges
    auto loads = path_table.getEdgeCongestion(graph);
    for (int e = 0; e < (int)loads.size(); ++e) {
        if (loads[e] > 0) {
            auto [u, v] = graph.getEdgeEndpoints(e);
            std::cout << "Edge (" << u << " -> " << v << "): load = " << loads[e] << "\n";
        }
    }

    // Get path length distribution
    auto dist = path_table.getPathLengthDistribution(source);
    std::cout << "Path lengths for source " << source << ":\n";
    for (const auto& [len, count] : dist) {
        std::cout << "  Length " << len << ": " << count << " paths\n";
    }

    // Print summary statistics
    path_table.printStats(graph);

    // Reconstruct linear flow table if needed
    LinearRoutingTable linear_table;
    path_table.reconstructLinearTable(linear_table, graph);
    linear_table.printFlows(graph);

    // Validate routing
    if (path_table.isValid(graph)) {
        std::cout << "Routing is valid!\n";
    }
}

*/

// Benefits of using PathRoutingTable:
//
// 1. CONGESTION ANALYSIS:
//    - Easily identify bottleneck edges
//    - Compute max congestion ratio
//    - Analyze load distribution
//
// 2. PATH ANALYSIS:
//    - Track actual routing paths used
//    - Measure path diversity
//    - Analyze path lengths
//    - Detect loops/cycles in paths
//
// 3. DEBUGGING:
//    - See exactly which edges carry flow
//    - Understand routing decisions
//    - Validate flow conservation
//
// 4. COMPARISON:
//    - Compare congestion of different routing strategies
//    - Analyze trade-offs (path length vs congestion)
//    - Benchmark algorithm performance
//
// 5. COMPATIBILITY:
//    - Can reconstruct to LinearRoutingTable for other analyses
//    - Works with all existing IGraph implementations
//    - Integrates with existing validation code

#endif

