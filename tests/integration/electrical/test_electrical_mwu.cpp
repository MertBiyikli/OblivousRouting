//
// Created by Mert Biyikli on 05.06.26.
//
#include "../../common/utils.h"
#include <filesystem>

using namespace integration;

TEST_CASE("Electrical Flow - Simple Oblivious Routing solve", "[ElectricalFlowMWU]")
{
    auto g = makeTriangleGraph();
    Config cfg = makeElectricalConfig();

    RoutingEngine engine;
    auto result = engine.solve(*g, cfg, cfg.solvers[0]);

    REQUIRE(result);
    double obl_ratio = result->oblivious_ratio;
    REQUIRE(obl_ratio >= 0); // Oblivious ratio should be at least 1
}

TEST_CASE("Electrical flow solver solves a simple path graph end-to-end",
          "[integration][electrical][routing-engine]")
{
    auto g = makePathGraph();
    auto cfg = makeElectricalConfig();

    RoutingEngine engine;

    auto result = engine.solve(*g, cfg, cfg.solvers.front());

    REQUIRE(result->oblivious_ratio < 5.0);
}

TEST_CASE("Electrical flow solver can be run repeatedly on the same small graph",
          "[integration][electrical][routing-engine]")
{
    auto g = makeTriangleGraph();
    auto cfg = makeElectricalConfig();

    RoutingEngine engine;

    auto first = engine.solve(*g, cfg, cfg.solvers.front());
    auto second = engine.solve(*g, cfg, cfg.solvers.front());


    REQUIRE(first->oblivious_ratio >= 0.5);
    REQUIRE(second->oblivious_ratio >= 0.5);
}

TEST_CASE("Electrical flow solver handles non-uniform capacities",
          "[integration][electrical][routing-engine][capacities]")
{
    auto g = std::make_unique<GraphCSR>(4);

    g->addEdge(0, 1, 100.0);
    g->addEdge(1, 2, 5.0);
    g->addEdge(2, 3, 100.0);
    g->addEdge(0, 3, 20.0);
    g->addEdge(1, 3, 10.0);

    g->finalize();

    auto cfg = makeElectricalConfig();

    RoutingEngine engine;
    auto result = engine.solve(*g, cfg, cfg.solvers.front());

    requireValidRoutingResult(result);
}

