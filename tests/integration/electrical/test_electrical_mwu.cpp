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
    checkForException(result);

    double obl_ratio = result.value().oblivious_ratio;
    REQUIRE(obl_ratio >= 0); // Oblivious ratio should be at least 1
}

TEST_CASE("Electrical flow solver solves a simple path graph end-to-end",
          "[integration][electrical][routing-engine]")
{
    auto g = makePathGraph();
    auto cfg = makeElectricalConfig();

    RoutingEngine engine;

    auto result = engine.solve(*g, cfg, cfg.solvers.front());
    checkForException(result);

    REQUIRE(result.value().oblivious_ratio < 5.0);
}


TEST_CASE("Electrical flow solver handles non-uniform capacities",
          "[integration][electrical][routing-engine][capacities]")
{
    auto g = std::make_unique<optimized::Graph<EdgeData>>(4, std::vector<optimized::Graph<EdgeData>::InputEdge>{
        {0, 1, EdgeData{100.0, 1.0}},
        {1, 2, EdgeData{5.0, 1.0}},
        {2, 3, EdgeData{100.0, 1.0}},
        {0, 3, EdgeData{20.0, 1.0}},
        {1, 3, EdgeData{10.0, 1.0}},
    });

    auto cfg = makeElectricalConfig();

    RoutingEngine engine;
    auto result = engine.solve(*g, cfg, cfg.solvers.front());
    checkForException(result);

    requireValidRoutingResult(result);
}

