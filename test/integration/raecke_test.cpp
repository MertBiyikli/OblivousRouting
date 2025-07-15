//
// Created by Mert Biyikli on 10.06.25.
//

#include <gtest/gtest.h>
#include "../../src/tree_based/raecke_tree_decomp.h"
#include "../../src/tree_based/raecke_transform.h"
#include "../../src/utils/hash.h"

// Test if graph is correctly assigned
TEST(RackeSolverTest, SetGraphCorrectlyStores) {
    auto graph = std::make_shared<Graph>(2);
    graph->addEdge(0, 1, 1.0);
    RackeObliviousRoutingSolver solver;
    solver.setGraph(graph);
    EXPECT_EQ(solver.getTrees().size(), 0); // Initially no trees
}

// Test createTreesAndLambda generates trees and lambdas
TEST(RackeSolverTest, CreateTreesAndLambdaGeneratesTrees) {
    auto graph = std::make_shared<Graph>(3);
    graph->addEdge(0, 1, 1.0);
    graph->addEdge(0, 2, 3.0);
    graph->addEdge(1, 2, 2.0);
    graph->createDistanceMatrix();

    RackeObliviousRoutingSolver solver;
    solver.setGraph(graph);
    solver.createTreesAndLambda();

    auto trees = solver.getTrees();
    auto lambdas = solver.getLambdas();
    EXPECT_GT(trees.size(), 0);
    EXPECT_GT(lambdas.size(), 0);

    double lambdaSum = 0.0;
    for (double l : lambdas) lambdaSum += l;
    EXPECT_NEAR(lambdaSum, 1.0, 1e-5);
}

// Test getTree runs MCCT and returns a tree
TEST(RackeSolverTest, GetTreeGeneratesValidTree) {
    auto graph = std::make_shared<Graph>(4);
    graph->addEdge(0, 1, 1.0);
    graph->addEdge(1, 2, 1.0);
    graph->addEdge(2, 3, 1.0);
    graph->createDistanceMatrix();

    RackeObliviousRoutingSolver solver;
    auto tree = solver.getTree(graph);

    EXPECT_NE(tree, nullptr);
    EXPECT_FALSE(tree->GetRoot()->GetVertices().empty());
}

// Test computeNewDistances returns positive distances
TEST(RackeSolverTest, ComputeNewDistancesBasicSanity) {
    auto graph = std::make_shared<Graph>(3);
    graph->addEdge(0, 1, 2.0);
    graph->addEdge(1, 2, 1.0);
    graph->createDistanceMatrix();

    RackeObliviousRoutingSolver solver;
    solver.setGraph(graph);
    solver.createTreesAndLambda();

    auto distances = solver.computeNewDistances(graph);
    EXPECT_GT(distances.size(), 0);
    for (const auto& [arc, dist] : distances) {
        EXPECT_GT(dist, 0.0);
    }
}


// Test rload calculation returns nonzero value
TEST(RackeSolverTest, RloadCalculationSanityCheck) {
    auto graph = std::make_shared<Graph>(3);
    graph->addEdge(0, 1, 1.0);
    graph->addEdge(0, 2, 1.0);
    graph->addEdge(1, 2, 1.0);
    graph->createDistanceMatrix();

    RackeObliviousRoutingSolver solver;
    solver.setGraph(graph);
    solver.createTreesAndLambda();

    double totalRload = solver.getRloadAllEdges(graph);
    EXPECT_GT(totalRload, 0.0);
}


TEST(RackeSolverTestExtended, LambdaWeightedConvexCombinationMatches) {
    auto graph = std::make_shared<Graph>(3);
    graph->addEdge(0, 1, 1.0);
    graph->addEdge(0, 2, 3.0);
    graph->addEdge(1, 2, 2.0);
    graph->createDistanceMatrix();

    RackeObliviousRoutingSolver solver;
    solver.setGraph(graph);
    solver.createTreesAndLambda();

    const auto& trees = solver.getTrees();
    const auto& lambdas = solver.getLambdas();

    ASSERT_EQ(trees.size(), lambdas.size());

    double weightedSum = 0.0;
    for (size_t i = 0; i < trees.size(); ++i) {
        int sumTree = 0;
        for (int l = 0; l < trees[i]->GetLayerCount(); ++l) {
            for (auto& node : trees[i]->GetLayerNodes(l)) {
                sumTree += node->GetVertices().size();
            }
        }
        weightedSum += lambdas[i] * sumTree;
    }

    EXPECT_GT(weightedSum, 0.0);  // There's some contribution from all trees
}


TEST(RackeSolverTestExtended, DistanceUpdateIsMonotonic) {
    auto graph = std::make_shared<Graph>(4);
    graph->addEdge(0, 1, 1.0);
    graph->addEdge(1, 2, 1.0);
    graph->addEdge(2, 3, 1.0);
    graph->createDistanceMatrix();

    RackeObliviousRoutingSolver solver;
    solver.setGraph(graph);

    std::unordered_map<std::pair<int, int>, double> beforeUpdate;
    for (int u : graph->GetVertices()) {
        for (const auto& arc : graph->neighbors(u)) {
            beforeUpdate[{u, arc->target}] = arc->distance;
        }
    }

    solver.createTreesAndLambda();

    for (int u : graph->GetVertices()) {
        for (const auto& arc : graph->neighbors(u)) {
            double oldDist = beforeUpdate[{u, arc->target}];
            double newDist = arc->distance;
            EXPECT_GE(newDist, oldDist);
        }
    }
}


TEST(RackeSolverStressTest, CycleGraph10Nodes) {
    const int N = 10;
    auto graph = std::make_shared<Graph>(N);
    for (int i = 0; i < N; ++i) {
        graph->addEdge(i, (i+1) % N, 1.0);
    }
    graph->createDistanceMatrix();

    RackeObliviousRoutingSolver solver;
    solver.setGraph(graph);
    solver.createTreesAndLambda();

    EXPECT_GT(solver.getTrees().size(), 0);
}


TEST(RackeSolverStressTest, StarGraphWithCenter) {
    auto graph = std::make_shared<Graph>(11);
    const int center = 0;
    for (int i = 1; i <= 10; ++i) {
        graph->addEdge(center, i, 1.0);
    }
    graph->createDistanceMatrix();

    RackeObliviousRoutingSolver solver;
    solver.setGraph(graph);
    solver.createTreesAndLambda();

    auto trees = solver.getTrees();
    EXPECT_GT(trees.size(), 0);
}


TEST(RackeSolverStressTest, BinaryTreeGraph) {
    int N = 15;
    auto graph = std::make_shared<Graph>(N);
    for (int i = 1; i < N; ++i) {
        int parent = (i - 1) / 2;
        graph->addEdge(parent, i, 1.0);
    }
    graph->createDistanceMatrix();

    RackeObliviousRoutingSolver solver;
    solver.setGraph(graph);
    solver.createTreesAndLambda();

    auto lambdas = solver.getLambdas();
    double total = std::accumulate(lambdas.begin(), lambdas.end(), 0.0);
    EXPECT_NEAR(total, 1.0, 1e-5);
}

TEST(RackeSolverIntegrationTest, FullPipelineSanityCheck) {
    // Step 1: Create a sample graph (triangle with a tail)
    auto graph = std::make_shared<Graph>(5);
    graph->addEdge(0, 1, 1.0);  // triangle part
    graph->addEdge(0, 2, 1.0);
    graph->addEdge(1, 2, 1.0);
    graph->addEdge(2, 3, 2.0);  // tail
    graph->addEdge(3, 4, 2.0);
    graph->createDistanceMatrix();

    // Step 2: Create the solver and assign the graph
    RackeObliviousRoutingSolver solver;
    solver.setGraph(graph);

    // Step 3: Run the full pipeline
    solver.createTreesAndLambda();

    // Step 4: Validate number of trees and lambda sum
    const auto& trees = solver.getTrees();
    const auto& lambdas = solver.getLambdas();

    ASSERT_FALSE(trees.empty());
    ASSERT_EQ(trees.size(), lambdas.size());

    double lambdaSum = std::accumulate(lambdas.begin(), lambdas.end(), 0.0);
    EXPECT_NEAR(lambdaSum, 1.0, 1e-6);

    // Step 5: Compute and normalize distances
    auto updatedDistances = solver.computeNewDistances(graph);
    solver.normalizeDistance(updatedDistances);

    // Step 6: Ensure all updated distances are >= 1 and finite
    for (const auto& [arc, dist] : updatedDistances) {
        EXPECT_TRUE(std::isfinite(dist));
        EXPECT_GE(dist, 1.0);
    }

    // Step 7: Ensure the r-load sum is meaningful
    double totalRload = solver.getRloadAllEdges(graph);
    EXPECT_GT(totalRload, 0.0);


    const auto& rloads = solver.GetEdge2Load().at(0);

    // Expect roughly 4 tree arcs (depends on MCCT structure)
    EXPECT_GE(rloads.size(), 3);
    EXPECT_LE(rloads.size(), 5);

    for (const auto& [arc, load] : rloads) {
        EXPECT_NE(arc, nullptr);
        EXPECT_GT(load, 0.0);
        EXPECT_LE(load, 1.0);  // only one original edge crosses each cut
        totalRload += load;
    }

    // Total rload across all tree arcs ≈ 4.0
    EXPECT_NEAR(totalRload, 4.0, 1e-2);
}


// Simple graph generator: ring topology for simplicity
std::shared_ptr<Graph> generateRingGraph(int n, double cap) {
    auto g = std::make_shared<Graph>(n);
    for (int i = 0; i < n; ++i) {
        int next = (i + 1) % n;
        g->addEdge(i, next, cap);
        g->addEdge(next, i, cap); // undirected
    }
    return g;
}

TEST(RaeckSolverIntegrationTest, RunSolve) {
    int n = 16; // number of nodes
    auto g = generateRingGraph(n, 10.0);

    // Set all edge distances to 1
    for (int u = 0; u < n; ++u)
        for (auto& e : g->neighbors(u))
            e->distance = 1.0;

    // Run Räcke's algorithm
    RackeObliviousRoutingSolver solver;
    solver.setGraph(g);
    solver.createTreesAndLambda();

    // Create oblivious routing solution
    RaeckeTransform transform;
    double lambda_sum = 0.0;
    DemandMap routing;

    for (size_t i = 0; i < solver.getTrees().size(); ++i) {
        double lambda_i = solver.getLambdas()[i];
        lambda_sum += lambda_i;
        double new_lambda = lambda_i / lambda_sum;
        routing = transform.addTree(solver.getTrees()[i], new_lambda, solver.getGraphs()[i]);
    }

    // Build all-pairs unit demands
    std::vector<std::pair<int, int>> demands;
    for (int u = 0; u < n; ++u)
        for (int v = 0; v < n; ++v)
            if (u != v) demands.emplace_back(u, v);

    // Compute arc loads
    std::unordered_map<std::pair<int, int>, double> arc_loads;
    for (const auto& [arc, dem_map] : routing) {
        for (const auto& [d, frac] : dem_map) {
            if (std::find(demands.begin(), demands.end(), d) != demands.end())
                arc_loads[arc] += frac;
        }
    }

    // Compute max congestion
    double max_cong = 0.0;
    for (const auto& [arc, load] : arc_loads) {
        auto& edges = g->neighbors(arc.first);
        for (const auto& edge : edges) {
            if (edge->target == arc.second) {
                max_cong = std::max(max_cong, load / edge->capacity);
                break;
            }
        }
    }

    EXPECT_LE(max_cong, std::log2(n));

    std::cout << "Max congestion on ring graph with n=" << n << " nodes: " << max_cong << "\n";
    std::cout << "Expected O(log n) bound: " << std::log2(n) << "\n";
}