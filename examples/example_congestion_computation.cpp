/**
 * Example: Using computeCongestion() to evaluate routing under demands
 *
 * This example demonstrates:
 * 1. computeCongestion(demand_map) - route specific demands through the routing table
 * 2. computeCongestionFromCapacities() - route edge capacities as demands
 * 3. Analysis of bottlenecks and load distribution
 */

#include "include/core/path_routing_table.h"
#include "include/data_structures/graph/Igraph.h"
#include <iostream>

// Example usage:
/*

void example_demand_based_congestion() {
    // Load graph
    IGraph graph;
    graph.loadFromLGF("network.lgf");

    // Create and populate path routing table
    PathRoutingTable path_table;
    path_table.init(graph);

    // Add routing paths (from some routing algorithm)
    // For example, from tree decomposition...
    // path_table.addPath(source, path_edges, flow);

    // =====================================================
    // METHOD 1: Compute congestion with custom demands
    // =====================================================

    // Create a demand map: source -> demand amount
    std::unordered_map<int, double> demands;
    demands[1] = 0.5;   // Source 1 sends 0.5 units of demand
    demands[2] = 0.3;   // Source 2 sends 0.3 units of demand
    demands[3] = 0.7;   // Source 3 sends 0.7 units of demand

    // Compute congestion under these demands
    double max_congestion = path_table.computeCongestion(graph, demands);
    std::cout << "Max congestion under custom demands: " << max_congestion << "\n";

    // Get detailed statistics
    path_table.printCongestionStats(graph, demands);

    // Get edge-by-edge loads
    auto loads = path_table.computeEdgeLoadsWithDemands(graph, demands);
    std::cout << "Edge loads:\n";
    for (int e = 0; e < (int)loads.size(); ++e) {
        if (loads[e] > 1e-6) {  // Only print non-zero
            auto [u, v] = graph.getEdgeEndpoints(e);
            std::cout << "  Edge (" << u << " -> " << v << "): " << loads[e] << "\n";
        }
    }

    // =====================================================
    // METHOD 2: Compute congestion from edge capacities
    // =====================================================

    // This creates a "worst case" scenario where each edge's capacity
    // is treated as a demand source from u to v
    double max_cong_from_capacities = path_table.computeCongestionFromCapacities(graph);
    std::cout << "\nMax congestion with capacity-based demands: "
              << max_cong_from_capacities << "\n";

    // Also print detailed stats for this scenario
    std::unordered_map<int, double> capacity_demands;
    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        auto [u, v] = graph.getEdgeEndpoints(e);
        double capacity = graph.getEdgeCapacity(e);
        capacity_demands[u] += capacity;
    }
    path_table.printCongestionStats(graph, capacity_demands);

    // =====================================================
    // METHOD 3: Compare different demand patterns
    // =====================================================

    std::cout << "\n========== DEMAND PATTERN COMPARISON ==========\n";

    // Pattern 1: Uniform demand
    std::unordered_map<int, double> uniform_demands;
    for (int src = 1; src < graph.getNumNodes(); ++src) {
        uniform_demands[src] = 1.0;  // Each source sends 1 unit
    }
    double cong_uniform = path_table.computeCongestion(graph, uniform_demands);
    std::cout << "Uniform demand congestion: " << cong_uniform << "\n";

    // Pattern 2: Concentrated demand
    std::unordered_map<int, double> concentrated_demands;
    concentrated_demands[1] = 5.0;  // Only source 1 sends large demand
    double cong_concentrated = path_table.computeCongestion(graph, concentrated_demands);
    std::cout << "Concentrated demand congestion: " << cong_concentrated << "\n";

    // Pattern 3: Skewed demand (heavy on early sources)
    std::unordered_map<int, double> skewed_demands;
    for (int src = 1; src < graph.getNumNodes(); ++src) {
        skewed_demands[src] = 1.0 / src;  // 1.0, 0.5, 0.33, ...
    }
    double cong_skewed = path_table.computeCongestion(graph, skewed_demands);
    std::cout << "Skewed demand congestion: " << cong_skewed << "\n";

    std::cout << "=============================================\n";
}

// Usage in algorithm evaluation:
void compare_routing_algorithms(const IGraph& graph) {
    // Algorithm A: Tree-based routing
    PathRoutingTable path_table_a;
    path_table_a.init(graph);
    // ... populate with tree-based paths ...

    // Algorithm B: Flow-based routing
    PathRoutingTable path_table_b;
    path_table_b.init(graph);
    // ... populate with flow-based paths ...

    // Create test demands
    std::unordered_map<int, double> demands;
    for (int src = 1; src < graph.getNumNodes(); ++src) {
        demands[src] = 1.0;
    }

    // Evaluate both algorithms
    double cong_a = path_table_a.computeCongestion(graph, demands);
    double cong_b = path_table_b.computeCongestion(graph, demands);

    std::cout << "Algorithm A congestion: " << cong_a << "\n";
    std::cout << "Algorithm B congestion: " << cong_b << "\n";

    if (cong_a < cong_b) {
        std::cout << "Algorithm A is better by "
                  << (cong_b - cong_a) / cong_b * 100 << "%\n";
    } else {
        std::cout << "Algorithm B is better by "
                  << (cong_a - cong_b) / cong_a * 100 << "%\n";
    }
}

// Example output:
/*

Max congestion under custom demands: 1.45

========== CONGESTION UNDER DEMAND ==========
Total sources with demand: 3
Total demand: 1.5

Edge Load Statistics:
  Max congestion: 1.45
  Min congestion (non-zero): 0.15
  Avg congestion (all edges): 0.18
  Std deviation: 0.32
  Used edges: 12 / 500

Top 10 most congested edges:
  1. Edge (3 -> 7): 1.45
  2. Edge (7 -> 0): 1.30
  3. Edge (2 -> 3): 1.15
  4. Edge (1 -> 2): 0.95
  5. Edge (5 -> 7): 0.85
  6. Edge (4 -> 5): 0.75
  7. Edge (6 -> 7): 0.65
  8. Edge (8 -> 9): 0.55
  9. Edge (9 -> 0): 0.45
  10. Edge (10 -> 11): 0.35
===========================================

Max congestion with capacity-based demands: 2.87

========== DEMAND PATTERN COMPARISON ==========
Uniform demand congestion: 3.42
Concentrated demand congestion: 5.15
Skewed demand congestion: 2.18
=============================================

*/

#endif

