/**
 * Integration Guide: Using PathRoutingTable with TreeTransform
 *
 * This shows how to track paths through tree decomposition
 * and evaluate congestion after cycle removal.
 */

#include "include/core/path_routing_table.h"
#include "include/core/routing_table.h"
#include "include/algorithms/mwu/oracle/tree/tree_transform.h"

// Example: Extending TreeTransform to track paths
/*

class PathTrackingTreeTransform : public TreeTransform {
private:
    PathRoutingTable path_table;
    std::vector<int> current_path;  // Track current path being built

public:
    // Modified applyFlow to also track paths
    void applyFlowWithTracking(
            const std::vector<int>& path,
            const std::vector<int>& childMembers,
            const std::unordered_set<int>& non_child_members,
            const std::unordered_set<int>& A,
            double lambda,
            LinearRoutingTable& table) {

        // For each source, track the path taken
        if (non_child_members.contains(root_node)) {
            for (int src : childMembers) {
                if (src == root_node) continue;

                // Current path is the sequence of edges
                // We can store this in path_table
                path_table.addPath(src, path, lambda);

                // Also add to linear table as before
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int e = graph.getEdgeId(path[i], path[i+1]);
                    table.addFlow(e, src, lambda);
                }
            }
        }

        // ... similar for destination case ...
    }

    // Get the path routing table after transformation
    PathRoutingTable& getPathTable() { return path_table; }
};

// Usage in algorithm:
void example_with_paths() {
    IGraph graph;
    graph.loadFromLGF("graph.lgf");

    PathTrackingTreeTransform tree_transform(graph);
    LinearRoutingTable routing_table;
    routing_table.init(graph);

    // Build HST and apply transformations
    // ... tree decomposition code ...

    // Get the paths that were taken
    PathRoutingTable paths = tree_transform.getPathTable();

    // Evaluate congestion BEFORE cycle removal
    double congestion_before = paths.evaluateCongestion(graph);
    std::cout << "Congestion before cycle removal: " << congestion_before << "\n";

    // Print congestion statistics
    auto stats = paths.getCongestionStats(graph);
    std::cout << "Congestion statistics:\n";
    std::cout << "  Max: " << stats["max"] << "\n";
    std::cout << "  Avg: " << stats["avg"] << "\n";
    std::cout << "  Std: " << stats["std_dev"] << "\n";

    // Perform cycle removal
    tree_transform.removeCyclesUsingTarjanSCCAlgorithm(routing_table);

    // Reconstruct paths after cycle removal
    PathRoutingTable paths_after;
    paths_after.init(graph);
    paths_after.reconstructLinearTable(routing_table, graph);

    // Evaluate congestion AFTER cycle removal
    double congestion_after = paths_after.evaluateCongestion(graph);
    std::cout << "Congestion after cycle removal: " << congestion_after << "\n";

    // Compare improvements
    std::cout << "Congestion reduction: "
              << (congestion_before - congestion_after) / congestion_before * 100
              << "%\n";
}

*/

// Key Integration Points:
//
// 1. PATH TRACKING IN TREE TRANSFORM
//    - Modify applyFlow() to record paths
//    - Store path sequence alongside flow amount
//
// 2. CONGESTION BEFORE/AFTER
//    - Measure congestion before cycle removal
//    - Measure congestion after cycle removal
//    - Analyze how cycle removal affects bottlenecks
//
// 3. IDENTIFY PROBLEMATIC PATHS
//    - Use edge loads to find heavy-traffic edges
//    - Correlate with cycle removal effectiveness
//
// 4. PERFORMANCE ANALYSIS
//    - Path length distribution shows routing complexity
//    - Compare algorithms by congestion profile
//
// 5. BOTTLENECK IDENTIFICATION
//    - Some edges may be heavily used after cycle removal
//    - Can guide further optimization

// Sample output when running with paths:
/*

========== PATH ROUTING TABLE STATISTICS ==========
Number of nodes: 100
Number of edges: 5000
Total paths: 198
Average path length: 4.5

Congestion Statistics:
  Max congestion: 12.45
  Min congestion (non-zero): 0.15
  Avg congestion (all edges): 0.39
  Std deviation: 1.23
  Used edges: 487 / 5000

Path length distribution (sample):
  Source 1:
    Length 3: 1 path(s)
    Length 5: 1 path(s)
===================================================

Congestion before cycle removal: 12.45
Congestion after cycle removal: 8.32
Congestion reduction: 33.1%

Bottleneck edges (top 5):
  Edge (15->23): load = 12.32
  Edge (42->18): load = 11.98
  Edge (67->89): load = 10.45
  Edge (34->56): load = 9.87
  Edge (78->22): load = 9.12

*/

#endif

