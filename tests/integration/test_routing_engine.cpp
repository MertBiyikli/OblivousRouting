#include "common/utils.h"
#include "core/routing_engine.h"

using namespace integration;

TEST_CASE("RoutingEngine solves small graph with FRT tree solver", "[integration][routing-engine][tree][frt]")
{
    runTreeSolverAndRequireValid(makeCycleGraph, SolverType::RAECKE_FRT_FLAT);
}

TEST_CASE("RoutingEngine solves small graph with CKR tree solver", "[integration][routing-engine][tree][ckr]")
{
    runTreeSolverAndRequireValid(makeCycleGraph, SolverType::RAECKE_CKR_FLAT);
}

TEST_CASE("RoutingEngine solves small graph with MST tree solver", "[integration][routing-engine][tree][mst]")
{
    runTreeSolverAndRequireValid(makeCycleGraph, SolverType::RAECKE_RANDOM_MST_FLAT);
}


TEST_CASE("RoutingEngine solves small graph with Electrical Flow solver",
          "[integration][routing-engine][electrical]")
{
    runElectricalSolverAndRequireValid(makeCycleGraph);
}

TEST_CASE("RoutingEngine solves small graph with Applegate-Cohen LP solver",
          "[integration][routing-engine][lp][applegate-cohen]")
{
    auto graph = makeTriangleGraph();

    Config cfg;
    cfg.solvers = {SolverType::LP_APPLEGATE_COHEN}; // adapt enum name
    cfg.graph_format = GraphFormat::CSR;

    RoutingEngine engine;
    auto result = engine.solve(*graph, cfg, cfg.solvers.front());

    requireValidRoutingResult(result);
}