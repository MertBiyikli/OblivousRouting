//
// Created by Mert Biyikli on 09.06.26.
//

#include "../../common/utils.h"
#include <catch2/catch_approx.hpp>
#include "algorithms/lp/lp_ac.h"
#include "algorithms/lp/lp_mcf.h"

#include "core/routing_table.h"

using namespace integration;

TEST_CASE("Applegate-Cohen LP solver runs on a triangle graph",
          "[integration][lp][applegate-cohen]")
{
    auto graph = makeTriangleGraph();

    LPSolver solver(*graph);

    AllPairRoutingTable table;
    table.init(*graph);
    solver.computeBasisFlows(table);

    //table.printFlows(*graph);
    requireValidAllPairRoutingTable(table, *graph);
}

TEST_CASE("Applegate-Cohen LP solver runs on a path graph",
          "[integration][lp][applegate-cohen]")
{
    auto graph = makePathGraph();

    LPSolver solver(*graph);

    AllPairRoutingTable table;
    table.init(*graph);
    solver.computeBasisFlows(table);

    requireValidAllPairRoutingTable(table, *graph);
}

TEST_CASE("CMMF solver runs on a triangle graph with one demand",
          "[integration][lp][mcf]")
{
    auto graph = makeTriangleGraph();

    CMMF_Solver solver(*graph);
    demands d;
    d.addDemand(0, 2, 1.0); // demand from node 0 to node 2 with value 1.0
    solver.AddDemandMap(d);
    auto offline_scheme = solver.solve();

    const double congestion = solver.getCongestionForPassedDemandMap();
    REQUIRE(std::isfinite(congestion));
    REQUIRE(congestion == Catch::Approx(0.5));
}