//
// Created by Mert Biyikli on 20.05.25.
//


#include "gtest/gtest.h"
#include "../../src/graph.h"
#include "../../src/tree_based/mcct/Node.h"
#include "../../src/tree_based/mcct/Tree.h"
#include "../../src/tree_based/mcct/mcct_derandomized_weighted_solver.h"
#include <stdexcept>
#include <algorithm>
#include <queue>
/*
// Unit Test for ContainsVertex function
TEST(NodeTest, ContainsVertexTest) {
    Node node;
    int v1(1);
    int v2(2);
    int v3(3);

    node.AddVertex(v1);
    node.AddVertex(v2);
    node.AddVertex(v3);

    // Test when vertex exists
    EXPECT_TRUE(node.ContainsVertex(1)); // Should return true
    EXPECT_TRUE(node.ContainsVertex(2)); // Should return true
    EXPECT_TRUE(node.ContainsVertex(3)); // Should return true

    // Test when vertex does not exist
    EXPECT_FALSE(node.ContainsVertex(4)); // Should return false
}



TEST(TreeTest, RootNodeCreation) {
    Tree tree;
    tree.CreateRoot(0);  // Create root node at layer 0

    // Test that the root node exists
    EXPECT_NE(tree.GetRoot(), nullptr);

    // Test the root node ID is set correctly
    EXPECT_EQ(tree.GetRoot()->GetId(), 0);  // The root should have ID 0 (the first node)
}


TEST(TreeTest, AddNodeToLayer) {
    Tree tree;
    tree.CreateRoot(0);  // Create root node at layer 0

    std::shared_ptr<Node> child1 = std::make_shared<Node>();
    child1->AddVertex(1);  // Add vertex to the child node

    tree.AddNodeToLayer(0, child1);  // Add the child node to layer 0

    // Test that the child node is added to layer 0
    auto layerNodes = tree.GetLayerNodes(0);
    EXPECT_EQ(layerNodes.size(), 2);  // There should be two nodes in layer 0 (root and child1)
    EXPECT_EQ(layerNodes.back()->GetId(), 1);  // The last added node should be child1 (ID 1)
}

TEST(TreeTest, CheckSingletonLayer) {
    Tree tree;

    std::vector<int> vertices = {0, 1, 2};

    tree.CreateRoot(0);  // Create root node at layer 0
    tree.GetRoot()->AddVertices(vertices);

    std::shared_ptr<Node> child1 = std::make_shared<Node>();
    child1->AddVertex(1);  // Add vertex to the child node
    tree.AddNodeToLayer(1, child1);  // Add child1 to layer 1

    // Test that the layer has only singleton nodes (nodes with exactly 1 vertex)
    EXPECT_TRUE(tree.IsLayerWithAllSingletons(1));  // Layer 1 should have all singleton nodes
}

TEST(TreeTest, GetLeafVertex) {
    Tree tree;
    tree.CreateRoot(0);  // Create root node at layer 0

    std::shared_ptr<Node> child1 = std::make_shared<Node>();
    child1->AddVertex(1);  // Add vertex 1 to child1
    tree.AddNodeToLayer(0, child1);  // Add child1 to layer 0

    std::shared_ptr<Node> leaf = tree.GetLeafVertex(1);  // Get leaf node with vertex ID 1
    EXPECT_NE(leaf, nullptr);  // Leaf should not be null
    EXPECT_EQ(leaf->GetId(), child1->GetId());  // The leaf should be the child1 node
}

TEST(TreeTest, CompareTrees) {
    Tree tree1;
    Tree tree2;
    tree1.SetId(1);
    tree2.SetId(2);

    // Test that the comparison operator works
    bool result = !(tree1 < tree2);
    EXPECT_TRUE(result);  // tree1 should be "less" than tree2 due to its smaller ID
}

TEST(MCCTSolverTest, AddDemandValidAndInvalid) {
    MCCTDerandomizedWeightedSolver solver;

    // Valid demand
    EXPECT_NO_THROW(solver.addDemand(1, 2, 5.0));

    // Invalid self-demand
    EXPECT_THROW(solver.addDemand(2, 2, 1.0), std::invalid_argument);

    // Invalid negative id
    EXPECT_THROW(solver.addDemand(-1, 2, 1.0), std::invalid_argument);
}

TEST(MCCTSolverTest, ComputeBetasBasic) {
    auto graph = std::make_shared<Graph>(3);
    graph->addEdge(0, 1, 4.0);
    graph->addEdge(0, 2, 8.0);
    graph->addEdge(1, 2, 6.0);

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(graph);
    graph->createDistanceMatrix();

    solver.computeBetas();
    auto betas = solver.getBetas();

    EXPECT_EQ(betas.size(), 0);
}

TEST(MCCTSolverTest, ResetRestoresCleanState) {
    MCCTDerandomizedWeightedSolver solver;
    solver.addDemand(0, 1, 3.5);
    solver.setBestB(2.0);

    solver.reset();

    EXPECT_TRUE(solver.getVerticesPermutation().empty());
    EXPECT_TRUE(solver.getBetas().empty());
    EXPECT_EQ(solver.getGraph(), nullptr);
}

TEST(MCCTSolverTest, ComputeExpectationReturnsFinite) {
    auto graph = std::make_shared<Graph>(3);
    graph->addEdge(0, 1, 2.0);
    graph->addEdge(0, 2, 4.0);
    graph->addEdge(1, 2, 1.0);
    graph->createDistanceMatrix();

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(graph);
    solver.addDemand(0, 1, 5.0);
    solver.addDemand(1, 2, 2.0);

    std::vector<int> perm = {0};
    std::set<int> unsettled = {1, 2};

    double result = solver.computeExpectation(2.0, perm, unsettled);
    EXPECT_TRUE(std::isfinite(result));
    //EXPECT_GT(result, 0);
}

TEST(MCCTSolverTest, BestTreeBuildsSomething) {
    auto graph = std::make_shared<Graph>(3);
    graph->addEdge(0, 1, 2.0);
    graph->addEdge(0, 2, 4.0);
    graph->addEdge(1, 2, 1.0);
    graph->createDistanceMatrix();

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(graph);
    solver.addDemand(0, 1, 5.0);
    solver.addDemand(1, 2, 2.0);

    Tree bestTree = solver.getBestTree();

    auto root = bestTree.GetRoot();
    EXPECT_NE(root, nullptr);
    EXPECT_FALSE(root->GetVertices().empty());

}

TEST(MCCTSolverTest, BestBetaAndPermutationShouldWork) {
    auto graph = std::make_shared<Graph>(3);
    graph->addEdge(0, 1, 1.0);
    graph->addEdge(0, 2, 2.0);
    graph->addEdge(1, 2, 1.0);
    graph->createDistanceMatrix();

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(graph);
    solver.addDemand(0, 2, 1.0);
    solver.addDemand(1, 2, 2.0);

    solver.findBestBetaAndPermutation();

    const auto& perm = solver.getVerticesPermutation();
    EXPECT_EQ(perm.size(), graph->GetVertices().size());
    EXPECT_GT(solver.getBestB(), 0.0);
}

TEST(MCCTSolverTest, ComputeBestTreeShouldGenerateLayers) {
    auto graph = std::make_shared<Graph>(3);
    graph->addEdge(0, 1, 1.0);
    graph->addEdge(0, 2, 2.0);
    graph->addEdge(1, 2, 1.5);
    graph->createDistanceMatrix();

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(graph);
    solver.addDemand(0, 2, 1.0);

    solver.findBestBetaAndPermutation();
    solver.computeBestTree();

    Tree tree = solver.getTree();
    EXPECT_GT(tree.GetLayerCount(), 0);
    EXPECT_NE(tree.GetRoot(), nullptr);
    EXPECT_FALSE(tree.GetRoot()->GetVertices().empty());
}

TEST(MCCTSolverTest, GetBestTreeIntegration) {
    auto graph = std::make_shared<Graph>(4);
    graph->addEdge(0, 1, 1.0);
    graph->addEdge(1, 2, 1.0);
    graph->addEdge(2, 3, 1.0);
    graph->createDistanceMatrix();

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(graph);
    solver.addDemand(0, 3, 5.0);
    solver.addDemand(1, 2, 3.0);

    Tree bestTree = solver.getBestTree();
    EXPECT_NE(bestTree.GetRoot(), nullptr);
    EXPECT_GT(bestTree.GetLayerCount(), 0);

    // Ensure there's at least some non-root node
    bool foundNonRoot = false;
    for (int i = 0; i < bestTree.GetLayerCount(); ++i) {
        if (!bestTree.GetLayerNodes(i).empty()) {
            for (const auto& n : bestTree.GetLayerNodes(i)) {
                if (n->GetId() != 0) foundNonRoot = true;
            }
        }
    }
    EXPECT_TRUE(foundNonRoot);
}


TEST(MCCTSolverUseCase, SmallGraphExample) {
    auto graph = std::make_shared<Graph>(4);
    graph->addEdge(0, 1, 1.0);
    graph->addEdge(1, 2, 2.0);
    graph->addEdge(2, 3, 3.0);
    graph->createDistanceMatrix();

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(graph);
    solver.addDemand(0, 2, 2.0);
    solver.addDemand(0, 3, 3.0);
    solver.addDemand(1, 3, 1.0);

    Tree tree = solver.getBestTree();

    // Check tree structure
    EXPECT_FALSE(tree.GetLayerNodes(0).empty());
    EXPECT_FALSE(tree.GetRoot()->GetVertices().empty());
}

TEST(MCCTDerandomizedWeightedSolverTest, AddDemandInvalidInputThrows) {
    MCCTDerandomizedWeightedSolver solver;
    EXPECT_THROW(solver.addDemand(-1, 2, 1.0), std::invalid_argument);
    EXPECT_THROW(solver.addDemand(2, 2, 1.0), std::invalid_argument);
}

TEST(MCCTDerandomizedWeightedSolverTest, TreeCoversAllVerticesExactlyOnce) {
    Graph g(3);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 1.0);
    g.addEdge(0, 2, 1.0);
    g.createDistanceMatrix();

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(std::make_shared<Graph>(g));
    solver.addDemand(0, 1, 1.0);
    solver.addDemand(1, 2, 1.0);
    solver.addDemand(2, 0, 1.0);

    for(int v : g.GetVertices()) {
        for(auto& u : g.neighbors(v)) {

            u->distance = 10;
        }
    }

    Tree tree = solver.getBestTree();

    std::set<int> seen;
    for (int layer = 0; layer <= tree.GetLayerCount(); ++layer) {
        for (const auto& node : tree.GetLayerNodes(layer)) {
            for (int v : node->GetVertices()) {
                ASSERT_FALSE(seen.count(v)) << "Vertex " << v << " occurs in multiple nodes.";
                seen.insert(v);
            }
        }
    }

    std::vector<int> vertices = g.GetVertices();
    EXPECT_EQ(seen.size(), vertices.size());
}

TEST(MCCTDerandomizedWeightedSolverTest, TreeHierarchyIsValid) {
    Graph g(4);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 1.0);
    g.addEdge(2, 3, 1.0);
    g.createDistanceMatrix();

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(std::make_shared<Graph>(g));
    solver.addDemand(0, 1, 1.0);
    solver.addDemand(1, 2, 1.0);
    solver.addDemand(2, 3, 1.0);

    Tree tree = solver.getBestTree();

    for (int layer = 1; layer <= tree.GetLayerCount(); ++layer) {
        for (auto& node : tree.GetLayerNodes(layer)) {
            ASSERT_NE(node->GetParent(), nullptr) << "Non-root node has no parent.";
            auto parent = node->GetParent();
            bool found = false;
            for (auto& child : parent->GetChildren()) {
                if (child == node) {
                    found = true;
                    break;
                }
            }
            EXPECT_TRUE(found) << "Parent does not list node as child.";
        }
    }
}

 */
/*
TEST(MCCTDerandomizedWeightedSolverTest, TreePreservesDistancesWithinLogStretch) {
    Graph g(4);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 1.0);
    g.addEdge(2, 3, 1.0);
    g.addEdge(0, 3, 3.0);
    g.createDistanceMatrix();

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(std::make_shared<Graph>(g));
    solver.addDemand(0, 3, 1.0);
    solver.addDemand(1, 2, 1.0);

    Tree tree = solver.getBestTree();
    auto arcMap = tree.GetTreeArcMap();

    // BFS to compute tree distances
    auto bfsTreeDistance = [&](int start) {
        std::map<int, double> dist;
        std::queue<std::pair<int, double>> q;
        std::set<int> visited;

        q.push({start, 0.0});
        while (!q.empty()) {
            auto [cur, d] = q.front(); q.pop();
            if (visited.count(cur)) continue;
            visited.insert(cur);
            dist[cur] = d;

            for (const auto& [edge, arc] : arcMap) {
                if (edge.first == cur && !visited.count(edge.second))
                    q.push({edge.second, d + arc.getLength()});
                if (edge.second == cur && !visited.count(edge.first))
                    q.push({edge.first, d + arc.getLength()});
            }
        }
        return dist;
    };

    auto vertices = g.GetVertices();
    double logn = std::log2(vertices.size());
    for (int u : vertices) {
        auto treeDistFromU = bfsTreeDistance(u);
        for (int v : vertices) {
            if (u == v) continue;
            double dG = g.GetDistanceMatrix().at(u).at(v);
            if (treeDistFromU.count(v)) {
                double dT = treeDistFromU[v];
                EXPECT_LE(dT, logn * dG * 5.0)
                                    << "Stretch too high from " << u << " to " << v << ": "
                                    << "graph=" << dG << " tree=" << dT;
            }
        }
    }
}
*/


class MCCTSolverTest{
public:
    static std::shared_ptr<Graph> makeTriangleGraph(int n = 3) {
        auto g = std::make_shared<Graph>(n);
        g->addEdge(0, 1, 1.0);
        g->addEdge(1, 2, 1.0);
        g->addEdge(0, 2, 2.0);
        for(int i = 0; i < 3; ++i) {
            for(auto& e: g->neighbors(i)) {
                e->distance = 10; // Set a constant distance for simplicity
            }
        }
        g->createDistanceMatrix();
        return g;
    }

    MCCTDerandomizedWeightedSolver solver;
};

TEST(MCCTSolverTest, ConstructorAndResetWorks) {
    MCCTDerandomizedWeightedSolver solver;
    solver.addDemand(0, 1, 1.0);
    solver.setBestB(2.0);
    solver.setVerticesPermutation({0, 1});
    solver.reset();

    EXPECT_TRUE(solver.getBetas().empty());
    EXPECT_TRUE(solver.getVerticesPermutation().empty());
    EXPECT_TRUE(solver.getGraph() == nullptr);
}

TEST(MCCTSolverTest, AddDemandRejectsInvalid) {
    MCCTDerandomizedWeightedSolver solver;
    EXPECT_THROW(solver.addDemand(-1, 2, 1.0), std::invalid_argument);
    EXPECT_THROW(solver.addDemand(1, 1, 1.0), std::invalid_argument);
}

TEST(MCCTSolverTest, ComputeBetasIncludesExpected) {
    MCCTDerandomizedWeightedSolver solver;
    auto g = MCCTSolverTest::makeTriangleGraph();
    solver.setGraph(g);
    solver.computeBetas();

    auto betas = solver.getBetas();
    EXPECT_FALSE(betas.empty());
    for (auto b : betas) {
        EXPECT_GT(b, 1.0);
    }
}

TEST(MCCTSolverTest, ComputeExpectationHandlesSimpleCase) {
    MCCTDerandomizedWeightedSolver solver;
    auto g = MCCTSolverTest::makeTriangleGraph();
    solver.setGraph(g);
    solver.addDemand(0, 1, 4.0);

    // Step 3: Set edge (0, 2) to have higher distance (20)
    int u = 0, v = 2;
    auto edge = std::find_if(
            g->neighbors(u).begin(), g->neighbors(u).end(),
            [v](const std::shared_ptr<Edge>& e) { return e->target == v; });
    if (edge != g->neighbors(u).end()) {
        (*edge)->distance = 20;
    }

    auto rev_edge = std::find_if(
            g->neighbors(v).begin(), g->neighbors(v).end(),
            [u](const std::shared_ptr<Edge>& e) { return e->target == u; });
    if (rev_edge != g->neighbors(v).end()) {
        (*rev_edge)->distance = 20;
    }

    std::vector<int> perm = {2};
    std::set<int> unsettled = {0, 1};
    g->createDistanceMatrix();
    solver.computeBetas();
    double beta = *solver.getBetas().begin();
    double cost = solver.computeExpectation(beta, perm, unsettled);

    EXPECT_GT(cost, 0.0);
}

TEST(MCCTSolverTest, FindBestBetaAndPermutationRunsCorrectly) {
    MCCTDerandomizedWeightedSolver solver;
    auto g = MCCTSolverTest::makeTriangleGraph();
    solver.setGraph(g);
    solver.addDemand(0, 1, 1.0);
    solver.addDemand(1, 2, 1.0);
    solver.findBestBetaAndPermutation();

    EXPECT_EQ(solver.getVerticesPermutation().size(), g->GetVertices().size());
    EXPECT_GT(solver.getBestB(), 0);
}


TEST(MCCTSolverTest, TreeHasSingletonLeaves) {
    MCCTDerandomizedWeightedSolver solver;
    auto g = MCCTSolverTest::makeTriangleGraph();
    solver.setGraph(g);
    solver.addDemand(0, 2, 1.0);
    solver.addDemand(1, 2, 1.0);
    Tree t = solver.getBestTree();

    bool hasSingleton = false;
    int SingletonCtr = 0;

    std::set<int> seen;
    std::queue<std::shared_ptr<Node>> Q;
    Q.push(t.GetRoot());
    seen.insert(t.GetRoot()->GetId());

    while (!Q.empty()) {
        auto node = Q.front();
        // Count leaves (nodes with a singleton vertex set)
        if (node->GetVertices().size() == 1) {
            hasSingleton = true;
            SingletonCtr++;
        }
        Q.pop();

        std::set<int> children_vertices;

        for (const auto& child : node->GetChildren()) {
            if (!child || seen.count(child->GetId())) continue;

            seen.insert(child->GetId());
            Q.push(child);
        }

    }
    EXPECT_TRUE(hasSingleton);
    EXPECT_EQ(SingletonCtr, g->numNodes());
}
TEST(MCCTSolverTest, TreeHasCorrectHierarchicalStructure) {
    // Step 1: Construct the graph
    auto g = std::make_shared<Graph>(4);
    g->addEdge(0, 1, 1.0);
    g->addEdge(1, 2, 1.0);
    g->addEdge(2, 3, 1.0);
    g->addEdge(0, 3, 1.0);

    // Step 2: Set all distances to 10
    for (int i = 0; i < 4; ++i) {
        for (auto& edge : g->neighbors(i)) {
            edge->distance = 10;
        }
    }

    // Step 3: Set edge (0, 3) to have higher distance (20)
    int u = 0, v = 3;
    auto edge = std::find_if(
            g->neighbors(u).begin(), g->neighbors(u).end(),
            [v](const std::shared_ptr<Edge>& e) { return e->target == v; });
    if (edge != g->neighbors(u).end()) {
        (*edge)->distance = 20;
    }

    auto rev_edge = std::find_if(
            g->neighbors(v).begin(), g->neighbors(v).end(),
            [u](const std::shared_ptr<Edge>& e) { return e->target == u; });
    if (rev_edge != g->neighbors(v).end()) {
        (*rev_edge)->distance = 20;
    }

    // Step 4: Run the solver
    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(g);
    solver.addDemand(0, 1, 1.0);
    solver.addDemand(1, 2, 1.0);

    Tree t = solver.getBestTree();

    // Step 5: Verify tree structure
    std::shared_ptr<Node> root = t.GetRoot();
    EXPECT_NE(root, nullptr);

    std::set<int> leaves_nodes;
    std::set<int> seen;
    std::queue<std::shared_ptr<Node>> Q;
    Q.push(root);
    seen.insert(root->GetId());

    while (!Q.empty()) {
        auto node = Q.front();
        Q.pop();

        std::set<int> children_vertices;

        for (const auto& child : node->GetChildren()) {
            if (!child || seen.count(child->GetId())) continue;

            seen.insert(child->GetId());
            Q.push(child);

            // Accumulate child vertices
            for(const auto& v : child->GetVertices()) {
                children_vertices.insert(v);
            }


            // Assert child vertices are subset of parent
            EXPECT_TRUE(std::includes(
                    node->GetVertices().begin(), node->GetVertices().end(),
                    child->GetVertices().begin(), child->GetVertices().end()));
        }

        if(!children_vertices.empty()) {
            // Assert parent and children collectively represent the same vertex set
            EXPECT_TRUE(std::includes(
                    node->GetVertices().begin(), node->GetVertices().end(),
                    children_vertices.begin(), children_vertices.end()));
            EXPECT_TRUE(std::includes(
                    children_vertices.begin(), children_vertices.end(),
                    node->GetVertices().begin(), node->GetVertices().end()));
        }
        // Count leaves (nodes with a singleton vertex set)
        if (node->GetVertices().size() == 1) {
            leaves_nodes.insert(*node->GetVertices().begin());
        }
    }

    // Final check: each original node should appear once as a leaf
    EXPECT_EQ(leaves_nodes.size(), g->numNodes());
}


TEST(MCCTSolverTest, TreeCostMatchesExpectation) {
    auto g = std::make_shared<Graph>(4);
    g->addEdge(0, 1, 1.0);
    g->addEdge(1, 2, 1.0);
    g->addEdge(2, 3, 1.0);
    g->addEdge(0, 3, 2.0);

    // Assign distances
    for (int i = 0; i < 4; ++i)
    for (auto& e : g->neighbors(i))
    e->distance = 10;

    for (auto& e : g->neighbors(0))
    if (e->target == 3) e->distance = 20;
    for (auto& e : g->neighbors(3))
    if (e->target == 0) e->distance = 20;

    MCCTDerandomizedWeightedSolver solver;
    solver.setGraph(g);
    solver.addDemand(0, 1, 1.0);
    solver.addDemand(1, 2, 1.0);

    // Store input before computing bestTree
    auto unsettled = std::set<int>(g->GetVertices().begin(), g->GetVertices().end());

    // Simulate computeExpectation
    solver.graph->createDistanceMatrix();  // Ensure distance matrix is built
    solver.computeBetas();

    double minCost = std::numeric_limits<double>::max();
    double bestB = -1;
    for (double beta : solver.getBetas()) {
        double cost = solver.computeExpectation(beta, solver.verticesPermutation, unsettled);
        if (cost < minCost) {
            minCost = cost;
            bestB = beta;
        }
    }

    // Set the selected bestB
    solver.setBestB(bestB);

    Tree tree = solver.getBestTree();

    // Re-traverse the tree and calculate expected cost from actual cuts
    double treeCost = 0.0;

    std::function<int(std::shared_ptr<Node>, int)> findLevel =
                [](std::shared_ptr<Node> node, int vertex) -> int {
                std::queue<std::pair<std::shared_ptr<Node>, int>> q;
                q.push({node, 0});
                while (!q.empty()) {
                    auto [curr, level] = q.front();
                    q.pop();
                    if (curr->ContainsVertex(vertex)) {
                        for (const auto& child : curr->GetChildren()) {
                            if (child->ContainsVertex(vertex)) {
                                q.push({child, level + 1});
                                break;
                            }
                        }
                        return level;
                    }
                }
                return -1;  // not found
            };

            for (const auto& [u, v2demand] : solver.idVertex2idVertex2demand) {
            for (const auto& [v, demand] : v2demand) {
            int lvl_u = findLevel(tree.GetRoot(), u);
            int lvl_v = findLevel(tree.GetRoot(), v);
            int lca_level = std::min(lvl_u, lvl_v);  // approx for now
            treeCost += std::pow(2, lca_level + 2) * demand;
        }
    }

    std::cout << "Expected Cost (computeExpectation) = " << minCost << "\n";
    std::cout << "Actual Tree Cost (from LCA levels) = " << treeCost << "\n";

    EXPECT_NEAR(minCost, treeCost, 1e-6);  // Allow small numerical error
}

