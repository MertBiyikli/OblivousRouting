# Test Build Configuration

This directory contains comprehensive unit tests for the ObliviousRouting project.

## Test Files

### Graph Data Structures
- **test_graph_csr.cpp**: Tests for Compressed Sparse Row (CSR) graph representation
  - Edge addition, finalization, and access
  - Shortest path algorithms
  - Bidirectional Dijkstra
  - Neighbor queries

- **test_graph_adj.cpp**: Tests for Adjacency List graph representation
  - Edge management and queries
  - Distance and capacity operations
  - Shortest path computations
  - Graph structure integrity

### Routing Algorithms
- **test_fast_ckr.cpp**: Tests for Fast CKR tree metric
  - Level partition computation
  - Cluster assignment
  - Radius settings
  - Large graph partitioning

- **test_frt.cpp**: Tests for FRT tree metric
  - Level partition with deterministic radii
  - Center assignment
  - Path-based clustering
  - Various delta values

### Network Flow
- **test_electrical_flow.cpp**: Tests for Electrical Flow (Laplacian Solver)
  - Laplacian matrix construction
  - Electrical potential computation
  - Different network topologies
  - Boundary conditions

### Routing Tables
- **test_routing_table.cpp**: Tests for routing table data structures
  - AllPairRoutingTable: Stores flows per edge and commodity pair
  - LinearRoutingTable: Stores flows per edge and source
  - Flow addition/removal
  - Multi-commodity routing

## Building Tests

### Prerequisites
- CMake 3.20+
- C++20 compiler
- Catch2 3.x (automatically fetched if not found)

### Build Instructions

#### Build with tests enabled:
```bash
cmake -DBUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release -S . -B build
cmake --build build -j$(nproc)
```

#### Run tests:
```bash
# Run all tests
cmake --build build --target run_tests

# Or directly with ctest
cd build
ctest --output-on-failure

# Run specific test category
ctest -R "GraphCSR" --output-on-failure

# Run with verbose output
ctest --output-on-failure -VV
```

## Test Coverage

The test suite covers:
- **Unit Tests**: Individual class functionality
- **Integration Tests**: Multi-class interactions
- **Edge Cases**: Empty inputs, large graphs, special values
- **Performance Tests**: Large graph handling

## Test Structure

Each test file follows Catch2 conventions:
```cpp
TEST_CASE("TestName", "[Category]") {
    // Setup
    
    // Test code
    REQUIRE(condition);
}
```

Tests are organized by category tags:
- `[GraphCSR]` - CSR graph tests
- `[GraphADJList]` - Adjacency list tests
- `[FastCKR]` - Fast CKR tests
- `[FRT]` - FRT tests
- `[LaplacianSolver]` - Electrical flow tests
- `[RoutingTable]` - Routing table tests

## Running Tests in CI/CD

Tests are automatically run in GitHub Actions:
1. Build with `-DBUILD_TESTS=ON`
2. Run `cmake --build build --target run_tests`
3. Test results are reported in the CI summary

## Continuous Integration

See `.github/workflows/ci.yml` for the complete CI pipeline that includes:
- Compilation checks
- Unit test execution
- Test result reporting

## Adding New Tests

To add new tests:

1. Create a new test file in `tests/` directory
2. Include necessary headers and Catch2
3. Write test cases using `TEST_CASE` macro
4. Add the file to `CMakeLists.txt` in the test executable
5. Run tests to verify

Example:
```cpp
#include <catch2/catch_test_macros.hpp>
#include "your_class.h"

TEST_CASE("YourClass - Feature", "[YourClass]") {
    YourClass obj(params);
    REQUIRE(obj.method() == expected_value);
}
```

## Debugging Tests

### Run a single test:
```bash
./build/unit_tests "[TestName]"
```

### Run tests matching a pattern:
```bash
./build/unit_tests "[Category]"
```

### Verbose output:
```bash
./build/unit_tests -v
```

### Generate XML report:
```bash
./build/unit_tests -r xml -o test_results.xml
```

## Test Metrics

Current test suite:
- **Total Tests**: 80+
- **Test Files**: 6
- **Coverage Areas**: Graphs, Algorithms, Routing, Network Flow

## Troubleshooting

### Tests fail to compile
- Ensure all headers are properly included
- Check that source files are added to CMakeLists.txt

### Tests fail at runtime
- Verify graph/algorithm parameters are valid
- Check memory allocations and bounds
- Look for floating point precision issues (use Approx)

### Catch2 not found
- CMake will automatically fetch Catch2 from GitHub
- Ensure internet connection is available during first build
- Or install Catch2 locally and set `Catch2_DIR`

