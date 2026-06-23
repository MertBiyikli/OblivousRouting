//
// Created by Mert Biyikli on 09.06.26.
//

#include <catch2/catch_test_macros.hpp>

#include <map>
#include <memory>
#include <vector>

#include "../../common/utils.h"

#include "algorithms/mwu/oracle/tree/frt/frt.h"
#include "algorithms/mwu/oracle/tree/fast_ckr/fast_ckr.h"
#include "algorithms/mwu/oracle/tree/tree_transform.h"

#include "../../../include/routing/routing_table.h"
#include "data_structures/hst/flat_hst.h"
#include "data_structures/hst/pointer_hst.h"


using namespace integration;

TEST_CASE("TreeTransform transforms a flat FRT tree into a linear routing table",
          "[integration][tree][transform][frt][flat-hst]")
{
    auto graph = makeCycleGraph();
    auto distances = unitDistances(*graph);

    FRT<FlatHST> oracle(*graph);
    auto hst = oracle.getTree(distances);

    LinearRoutingTable table;
    table.init(*graph);

    TreeTransform transform(*graph);
    std::map<std::pair<int, int>, CachedPath> path_cache;

    TreeIteration<FlatHST> iteration(std::move(hst), distances, 1.0);
    transform.transform(iteration, table, path_cache);

    requireNonEmptyLinearTable(table, *graph);
}

TEST_CASE("TreeTransform transforms a pointer FRT tree into a linear routing table",
          "[integration][tree][transform][frt][pointer-hst]")
{
    auto graph = makeCycleGraph();
    auto distances = unitDistances(*graph);

    FRT<std::shared_ptr<HSTNode>> oracle(*graph);
    auto hst = oracle.getTree(distances);

    LinearRoutingTable table;
    table.init(*graph);

    TreeTransform transform(*graph);
    std::map<std::pair<int, int>, CachedPath> path_cache;

    TreeIteration<std::shared_ptr<HSTNode>> iteration(std::move(hst), distances, 1.0);
    transform.transform(iteration, table, path_cache);

    requireNonEmptyLinearTable(table, *graph);
}