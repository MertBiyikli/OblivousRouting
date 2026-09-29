//
// Created by Mert Biyikli on 09.06.26.
//

#ifndef OBLIVIOUSROUTING_UTILS_H
#define OBLIVIOUSROUTING_UTILS_H

#pragma once

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <filesystem>
#include <memory>
#include <string>

#include "../../include/core/errors.h"
#include "../../include/routing/routing_engine.h"
#include "../../include/routing/routing_result.h"
#include "../../include/data_structures/graph/graph.h"
#include "../../include/io/parse_argument_io.h"

namespace integration {
    inline std::filesystem::path projectSourceDir()
    {
        return std::filesystem::path(PROJECT_SOURCE_DIR);
    }

    inline std::filesystem::path projectBinaryDir()
    {
        return std::filesystem::path(PROJECT_BINARY_DIR);
    }

    inline std::filesystem::path testDataRoot()
    {
        return projectSourceDir() / "tests" / "tiny_dataset";
    }

    inline std::filesystem::path tinyLgfDataset()
    {
        return testDataRoot() / "data" / "tiny_1221.lgf";
    }

    inline std::filesystem::path Backbone_1239_LgfDataset()
    {
        return testDataRoot() / "data" / "tiny_1239.lgf";
    }

    inline void requireExistingFile(const std::filesystem::path& path)
    {
        INFO("PROJECT_SOURCE_DIR = " << PROJECT_SOURCE_DIR);
        INFO("PROJECT_BINARY_DIR = " << PROJECT_BINARY_DIR);
        INFO("Current path        = " << std::filesystem::current_path().string());
        INFO("Requested file      = " << path.string());

        REQUIRE(std::filesystem::exists(path));
        REQUIRE(std::filesystem::is_regular_file(path));
    }

    inline std::unique_ptr<optimized::Graph<EdgeData>> makePathGraph()
    {
        std::vector<optimized::Graph<EdgeData>::InputEdge> edges = {
            {0, 1, EdgeData{10.0, 1.0}},
            {1, 2, EdgeData{10.0, 1.0}},
            {2, 3, EdgeData{10.0, 1.0}},
        };
        return std::make_unique<optimized::Graph<EdgeData>>(4, edges);
    }

    inline std::unique_ptr<optimized::Graph<EdgeData>> makeTriangleGraph()
    {
        std::vector<optimized::Graph<EdgeData>::InputEdge> edges = {
            {0, 1, EdgeData{1.0, 1.0}},
            {1, 2, EdgeData{1.0, 1.0}},
            {0, 2, EdgeData{1.0, 1.0}},
        };
        return std::make_unique<optimized::Graph<EdgeData>>(3, edges);
    }

    inline std::unique_ptr<optimized::Graph<EdgeData>> makeCycleGraph()
    {
        std::vector<optimized::Graph<EdgeData>::InputEdge> edges = {
            {0, 1, EdgeData{10.0, 1.0}},
            {1, 2, EdgeData{10.0, 1.0}},
            {2, 3, EdgeData{10.0, 1.0}},
            {3, 0, EdgeData{10.0, 1.0}},
        };
        return std::make_unique<optimized::Graph<EdgeData>>(4, edges);
    }

    inline std::unique_ptr<optimized::Graph<EdgeData>> makeCompleteGraphK4()
    {
        std::vector<optimized::Graph<EdgeData>::InputEdge> edges = {
            {0, 1, EdgeData{10.0, 1.0}},
            {0, 2, EdgeData{10.0, 1.0}},
            {0, 3, EdgeData{10.0, 1.0}},
            {1, 2, EdgeData{10.0, 1.0}},
            {1, 3, EdgeData{10.0, 1.0}},
            {2, 3, EdgeData{10.0, 1.0}},
        };
        return std::make_unique<optimized::Graph<EdgeData>>(4, edges);
    }

    inline std::unique_ptr<optimized::Graph<EdgeData>> makeNonUniformCapacityGraph()
    {
        std::vector<optimized::Graph<EdgeData>::InputEdge> edges = {
            {0, 1, EdgeData{100.0, 1.0}},
            {1, 2, EdgeData{5.0, 1.0}},
            {2, 3, EdgeData{100.0, 1.0}},
            {0, 3, EdgeData{20.0, 1.0}},
            {1, 3, EdgeData{10.0, 1.0}},
        };
        return std::make_unique<optimized::Graph<EdgeData>>(4, edges);
    }

    template <typename ResultPtr>
    inline void requireValidRoutingResult(const ResultPtr& result)
    {
        REQUIRE(std::isfinite(result.value().oblivious_ratio));
        REQUIRE(result.value().status == ResultStatus::OK);
    }

    template <typename T>
    inline void checkForException(const Result<T>& result) {
        if (!result) {
            FAIL("Failed running test: " << result.error().message);
        }
    }

    inline Config makeElectricalConfig()
    {
        Config cfg;
        cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
        cfg.graph_format = GraphFormat::CSR;
        cfg.evaluate_demand_models = false;

        return cfg;
    }

    inline Config makeTreeConfig(SolverType solver) {
        Config cfg;
        cfg.solvers = {solver};
        cfg.graph_format = GraphFormat::CSR;
        return cfg;
    }

    template <typename GraphFactory>
    inline void runTreeSolverAndRequireValid(GraphFactory&& graph_factory, SolverType solver_type) {
        auto graph = graph_factory();
        auto cfg = makeTreeConfig(solver_type);

        RoutingEngine engine;
        auto result = engine.solve(*graph, cfg, cfg.solvers.front());

        requireValidRoutingResult(result);
    }

    template <typename GraphFactory>
    inline void runElectricalSolverAndRequireValid(GraphFactory&& graph_factory) {
        auto graph = graph_factory();
        auto cfg = makeElectricalConfig();

        RoutingEngine engine;
        auto result = engine.solve(*graph, cfg, cfg.solvers.front());

        requireValidRoutingResult(result);
    }

    inline std::pair<int, int> normalizedEdge(int u, int v) {
        if (u > v) std::swap(u, v);
        return {u, v};
    }


    inline std::vector<double> unitDistances(const optimized::Graph<EdgeData>& graph)
    {
        return std::vector<double>(graph.getNumDirectedEdges(), 1.0);
    }

    inline void requireValidFlatHST(const FlatHST& tree, int expected_vertices)
    {
        REQUIRE(tree.size() >= expected_vertices);

        const auto root_members = tree.memberRange(tree.root());
        REQUIRE(static_cast<int>(root_members.size()) == expected_vertices);
    }

    inline void requireValidPointerHST(const std::shared_ptr<HSTNode>& root, int expected_vertices)
    {
        REQUIRE(root != nullptr);
        REQUIRE(!root->getChildren().empty());
        REQUIRE(static_cast<int>(root->getMembers().size()) == expected_vertices);
    }

    inline void requireNonEmptyLinearTable(const LinearRoutingTable& table, const optimized::Graph<EdgeData>& graph)
    {
        REQUIRE(table.getSize() != 0);
        REQUIRE(table.getNumNodes() == graph.getNumNodes());
    }


    inline void requireValidAllPairRoutingTable(const AllPairRoutingTable& table, const optimized::Graph<EdgeData>& graph)
    {
        REQUIRE(table.getSize() != 0);
        REQUIRE(table.getNumNodes() == graph.getNumNodes());
        REQUIRE(table.isValid(graph));
    }

};
#endif //OBLIVIOUSROUTING_UTILS_H