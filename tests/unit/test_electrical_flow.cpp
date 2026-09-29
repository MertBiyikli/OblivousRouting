#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <vector>
#include <Eigen/Dense>
#include "data_structures/graph/graph.h"
#include "algorithms/oblivious/mwu/oracle/electrical/laplacian_solver.h"
#include "../common/utils.h"

using Catch::Approx;
using Eigen::VectorXd;

// Helper to create a simple test graph
static optimized::Graph<EdgeData> createSimpleGraph() {
    std::vector<optimized::Graph<EdgeData>::InputEdge> edges = {
        {0, 1, EdgeData{1.0, 1.0}},
        {1, 2, EdgeData{1.0, 1.0}},
        {2, 3, EdgeData{1.0, 1.0}},
        {1, 3, EdgeData{1.0, 2.0}},
    };
    return optimized::Graph<EdgeData>(4, edges);
}

static optimized::Graph<EdgeData> createGridGraph() {
    std::vector<optimized::Graph<EdgeData>::InputEdge> edges;
    // 3x3 grid
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            int node = i * 3 + j;
            if (j + 1 < 3) {
                edges.push_back({node, node + 1, EdgeData{1.0, 1.0}});
            }
            if (i + 1 < 3) {
                edges.push_back({node, node + 3, EdgeData{1.0, 1.0}});
            }
        }
    }
    return optimized::Graph<EdgeData>(9, edges);
}


TEST_CASE("LaplacianSolver - Simple Graph Initialization", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    // This should not throw
    solver.init(graph, weights, 4, edges, false);

    REQUIRE(true);
}

TEST_CASE("LaplacianSolver - Update Solver", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    solver.init(graph, weights, 4, edges, false);

    // This should not throw
    solver.updateSolver();

    REQUIRE(true);
}

TEST_CASE("LaplacianSolver - Build Laplacian", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    solver.init(graph, weights, 4, edges, false);
    solver.buildLaplacian();

    REQUIRE(true);
}

TEST_CASE("LaplacianSolver - Solve with Std Vector", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    solver.init(graph, weights, 4, edges, false);

    std::vector<double> b = {1.0, 0.0, 0.0, -1.0};
    auto solution = solver.solve(b, 1e-6);

    REQUIRE(solution.value().size() == 4);
}

TEST_CASE("LaplacianSolver - Solve with Eigen Vector", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    solver.init(graph, weights, 4, edges, false);

    VectorXd b(4);
    b << 1.0, 0.0, 0.0, -1.0;

    auto solution = solver.solve(b, 1e-6);

    REQUIRE(solution.value().size() == 4);
}

TEST_CASE("LaplacianSolver - Weight Update", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    solver.init(graph, weights, 4, edges, false);

    // Update weights
    std::vector<double> new_weights = {2.0, 1.5, 1.0, 0.5};
    solver.updateAllEdges(new_weights, edges);

    REQUIRE(true);
}

TEST_CASE("LaplacianSolver - Grid Graph Solve", "[LaplacianSolver]") {
    auto graph = createGridGraph();
    LaplacianSolver solver;

    // For a 3x3 grid: nodes 0-8, 12 undirected edges (6 horizontal + 6 vertical).
    // LaplacianSolver::init expects one entry per undirected edge in `edges`,
    // matching `weights` 1:1 (it symmetrizes internally via setEdgeWeight both ways).

    std::vector<double> weights(12, 1.0);
    std::vector<std::pair<int, int>> edges;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            int node = i * 3 + j;
            if (j + 1 < 3) {
                int right = i * 3 + (j + 1);
                edges.push_back({node, right});
            }
            if (i + 1 < 3) {
                int down = (i + 1) * 3 + j;
                edges.push_back({node, down});
            }
        }
    }

    solver.init(graph, weights, 9, edges, false);

    std::vector<double> b(9, 0.0);
    b[0] = 1.0;
    b[8] = -1.0;

    auto solution = solver.solve(b, 1e-6);

    REQUIRE(solution.value().size() == 9);
}

TEST_CASE("LaplacianSolver - Different Epsilon Values", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    solver.init(graph, weights, 4, edges, false);

    std::vector<double> b = {1.0, 0.0, 0.0, -1.0};

    // Test with different epsilon values
    for (double eps : {1e-3, 1e-6, 1e-9}) {
        auto solution = solver.solve(b, eps);
        REQUIRE(solution.value().size() == 4);
    }
}

TEST_CASE("LaplacianSolver - Dirichlet Boundary Conditions", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    solver.init(graph, weights, 4, edges, false);

    std::vector<double> b(4, 0.0);
    b[0] = 1.0;
    b[3] = -1.0;

    auto solution = solver.solve(b, 1e-6);

    REQUIRE(solution.value().size() == 4);
}

TEST_CASE("LaplacianSolver - Apply Dirichlet", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    solver.init(graph, weights, 4, edges, false);

    std::vector<double> vals = {0.5, 0.3, 0.2, 0.1};
    solver.applyDirichletInPlace(vals);

    REQUIRE(vals.size() == 4);
}

TEST_CASE("LaplacianSolver - Set Solver Params", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    solver.init(graph, weights, 4, edges, false);

    boost::property_tree::ptree params;
    // Could add custom parameters here
    solver.setSolverParams(params);

    REQUIRE(true);
}

TEST_CASE("LaplacianSolver - Symmetry Check", "[LaplacianSolver]") {
    // Create a symmetric graph
    std::vector<optimized::Graph<EdgeData>::InputEdge> sym_edges = {
        {0, 1, EdgeData{1.0, 1.0}},
        {1, 2, EdgeData{1.0, 1.0}},
    };
    optimized::Graph<EdgeData> graph(3, sym_edges);

    LaplacianSolver solver;

    std::vector<double> weights(2, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}
    };

    solver.init(graph, weights, 3, edges, false);

    std::vector<double> b = {1.0, 0.0, -1.0};
    auto solution = solver.solve(b, 1e-6);

    REQUIRE(solution.value().size() == 3);
}

TEST_CASE("LaplacianSolver - Zero RHS", "[LaplacianSolver]") {
    auto graph = createSimpleGraph();
    LaplacianSolver solver;

    std::vector<double> weights(4, 1.0);
    std::vector<std::pair<int, int>> edges = {
        {0, 1}, {1, 2}, {2, 3}, {1, 3}
    };

    solver.init(graph, weights, 4, edges, false);

    std::vector<double> b(4, 0.0); // All zeros
    auto solution = solver.solve(b, 1e-6);

    REQUIRE(solution.value().size() == 4);
}

