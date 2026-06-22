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

#include "../../include/core/utils.h"
#include "../../include/core/routing_engine.h"
#include "../../include/data_structures/graph/graph_csr.h"
#include "../../include/io/parse_argurment_io.h"

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

    inline std::unique_ptr<GraphCSR> makePathGraph()
    {
        auto g = std::make_unique<GraphCSR>(4);

        g->addEdge(0, 1, 10.0);
        g->addEdge(1, 2, 10.0);
        g->addEdge(2, 3, 10.0);

        g->finalize();
        return g;
    }

    inline std::unique_ptr<GraphCSR> makeTriangleGraph()
    {
        auto g = std::make_unique<GraphCSR>(3);

        g->addEdge(0, 1, 1.0);
        g->addEdge(1, 2, 1.0);
        g->addEdge(0, 2, 1.0);

        g->finalize();
        return g;
    }

    inline std::unique_ptr<GraphCSR> makeCycleGraph()
    {
        auto g = std::make_unique<GraphCSR>(4);

        g->addEdge(0, 1, 10.0);
        g->addEdge(1, 2, 10.0);
        g->addEdge(2, 3, 10.0);
        g->addEdge(3, 0, 10.0);

        g->finalize();
        return g;
    }

    inline std::unique_ptr<GraphCSR> makeCompleteGraphK4()
    {
        auto g = std::make_unique<GraphCSR>(4);

        g->addEdge(0, 1, 10.0);
        g->addEdge(0, 2, 10.0);
        g->addEdge(0, 3, 10.0);
        g->addEdge(1, 2, 10.0);
        g->addEdge(1, 3, 10.0);
        g->addEdge(2, 3, 10.0);

        g->finalize();
        return g;
    }

    inline std::unique_ptr<GraphCSR> makeNonUniformCapacityGraph()
    {
        auto g = std::make_unique<GraphCSR>(4);

        g->addEdge(0, 1, 100.0);
        g->addEdge(1, 2, 5.0);
        g->addEdge(2, 3, 100.0);
        g->addEdge(0, 3, 20.0);
        g->addEdge(1, 3, 10.0);

        g->finalize();
        return g;
    }

    template <typename ResultPtr>
    inline void requireValidRoutingResult(const ResultPtr& result)
    {
        REQUIRE(result);
        REQUIRE(std::isfinite(result->oblivious_ratio));
        REQUIRE(result->status == ResultStatus::OK);
    }

    inline Config makeElectricalConfig()
    {
        Config cfg;
        cfg.solvers = {SolverType::ELECTRICAL_SKETCHING};
        cfg.graph_format = GraphFormat::CSR;

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


    inline std::vector<double> unitDistances(const IGraph& graph)
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

    inline void requireNonEmptyLinearTable(const LinearRoutingTable& table, const IGraph& graph)
    {
        REQUIRE(table.getSize() != 0);
        REQUIRE(table.getNumNodes() == graph.getNumNodes());
    }


    inline void requireValidAllPairRoutingTable(const AllPairRoutingTable& table, const IGraph& graph)
    {
        REQUIRE(table.getSize() != 0);
        REQUIRE(table.getNumNodes() == graph.getNumNodes());
        REQUIRE(table.isValid(graph));
    }

};
#endif //OBLIVIOUSROUTING_UTILS_H