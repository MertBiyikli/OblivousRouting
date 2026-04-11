# PathRoutingTable: Path-based Routing with Congestion Analysis

## Overview

`PathRoutingTable` is a new routing table data structure that stores **actual routing paths** instead of just edge flows. This enables:

1. **Congestion Evaluation** - Measure load on edges and network bottlenecks
2. **Path Analysis** - Track path diversity, lengths, and characteristics
3. **Flow Reconstruction** - Convert paths back to edge-flow representation
4. **Validation** - Verify flow conservation and routing correctness

## Key Features

### 1. Path Storage
Each routing stores paths as sequences of edge IDs with associated flow amounts:

```cpp
struct PathFlow {
    std::vector<int> path;  // Sequence of edge IDs
    double flow;             // Flow amount on this path
};

// For each source: list of paths
std::unordered_map<int, std::vector<PathFlow>> source_paths;
```

### 2. Congestion Evaluation

**Max Congestion:**
```cpp
double max_congestion = path_table.evaluateCongestion(graph);
```
Returns the maximum load on any edge in the network.

**Edge-by-Edge Loads:**
```cpp
std::vector<double> loads = path_table.getEdgeCongestion(graph);
```
Returns load on each edge, useful for identifying bottlenecks.

**Detailed Statistics:**
```cpp
auto stats = path_table.getCongestionStats(graph);
// Returns: max, min, avg, std_dev, used_edges, total_edges
```

### 3. Path Analysis

**Path Distribution:**
```cpp
auto dist = path_table.getPathLengthDistribution(source);
// Returns map from path_length -> count
```

**Global Metrics:**
```cpp
int total_paths = path_table.getTotalPathCount();
double avg_length = path_table.getAveragePathLength();
```

### 4. Flow Reconstruction

Convert paths back to edge-flow representation:

```cpp
LinearRoutingTable linear_table;
path_table.reconstructLinearTable(linear_table, graph);
```

This allows integration with existing flow-based analyses.

## Usage Example

```cpp
// Create and initialize
PathRoutingTable path_table;
path_table.init(graph);

// Add routing paths
int source = 1;
std::vector<int> path = {e0, e1, e2};  // Edge sequence
path_table.addPath(source, path, 0.5);  // Flow amount

// Evaluate congestion
double max_cong = path_table.evaluateCongestion(graph);
auto stats = path_table.getCongestionStats(graph);

std::cout << "Max congestion: " << stats["max"] << "\n";
std::cout << "Avg congestion: " << stats["avg"] << "\n";
std::cout << "Used edges: " << stats["used_edges"] << "\n";

// Print detailed stats
path_table.printStats(graph);
```

## Congestion Metrics Explained

### Maximum Congestion
The load on the most congested edge - useful for understanding the bottleneck:
```
max_congestion = max(flow on edge e) for all edges e
```

### Average Congestion
Average load across all edges (including unused):
```
avg = Σ(flow on edge e) / number_of_edges
```

### Std Deviation
Measure of load variance - lower means more uniform distribution:
```
std_dev = sqrt(Σ(load - avg)²) / num_edges
```

### Edge Utilization
Percentage of edges with non-zero load:
```
utilization = used_edges / total_edges
```

## Data Structure Files

### Header File
- **Location:** `include/core/path_routing_table.h`
- **Size:** ~150 lines
- **Classes:** `PathRoutingTable`

### Implementation
- **Location:** `source/core/path_routing_table.cpp`
- **Size:** ~280 lines
- **Methods:** 12+ public methods

### Example Usage
- **Location:** `examples/example_path_routing.cpp`
- **Demonstrates:** Adding paths, evaluating congestion, analyzing statistics

## Integration with Existing Code

### With LinearRoutingTable
```cpp
PathRoutingTable path_table;
// ... populate with paths ...

LinearRoutingTable linear_table;
path_table.reconstructLinearTable(linear_table, graph);
// Now use linear_table with existing code
```

### With IGraph
Works with any `IGraph` implementation:
- Standard graphs
- Weighted graphs
- Directed/undirected graphs

### With RoutingTable Interface
Inherits from `RoutingTable` base class:
- `init(const IGraph& g)`
- `isValid(const IGraph& g)`
- `printFlows(const IGraph& g)`

## Performance Characteristics

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| Add path | O(path_length) | Stores path sequence |
| Evaluate congestion | O(total_paths × avg_length) | Sums flows per edge |
| Get edge loads | O(total_paths × avg_length) | Pre-computes all edge loads |
| Reconstruct linear | O(total_paths × avg_length) | One-time conversion |
| Path statistics | O(total_paths) | Efficient aggregation |

## Memory Usage

For a network with:
- N nodes
- M edges  
- P total paths (across all sources)
- L average path length

**Memory:** O(P × L + M) for storing paths and computed statistics

## Use Cases

### 1. Algorithm Evaluation
Compare different routing algorithms by their congestion profiles:
```cpp
// Algorithm A
path_table_a.evaluateCongestion(graph);
// Algorithm B
path_table_b.evaluateCongestion(graph);
// Compare results
```

### 2. Bottleneck Analysis
Identify most congested edges:
```cpp
auto loads = path_table.getEdgeCongestion(graph);
// Find edges with load > threshold
```

### 3. Path Diversity
Analyze how traffic is distributed across paths:
```cpp
auto dist = path_table.getPathLengthDistribution(source);
// Longer paths = more diverse routing
```

### 4. Flow Verification
Ensure routing is valid before running algorithms:
```cpp
if (!path_table.isValid(graph)) {
    std::cerr << "Invalid routing!\n";
}
```

## Future Extensions

Possible enhancements:

1. **Per-Edge Statistics**
   - Track which sources use each edge
   - Compute edge utilization ratio

2. **Path Constraints**
   - Maximum path length enforcement
   - Disjoint path families

3. **Dynamic Updates**
   - Add/remove paths efficiently
   - Incremental congestion updates

4. **Visualization**
   - Export to format for network visualization
   - Heatmaps of edge loads

5. **Multi-Commodity**
   - Support all-pairs routing
   - Track flows between any source-destination pair

## Notes

- Paths are stored as **directed edge sequences**
- Flow conservation is validated per source
- Edge loads are computed on-demand from paths
- Zero-load edges are tracked but ignored in min/avg calculations
- Congestion is total flow on an edge (sum of all commodity flows)

---

**Created:** April 11, 2026  
**Author:** Mert Biyikli  
**Status:** Complete and tested

