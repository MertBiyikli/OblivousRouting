//
// Created by Mert Biyikli on 09.06.26.
//
#include <catch2/catch_test_macros.hpp>

#include <memory>
#include <vector>

#include "../common/utils.h"

#include "algorithms/mwu/oracle/tree/frt/frt.h"
#include "algorithms/mwu/oracle/tree/fast_ckr/fast_ckr.h"
#include "algorithms/mwu/oracle/tree/mst/mst_oracle.h"

#include "data_structures/hst/pointer_hst.h"
#include "data_structures/hst/flat_hst.h"

using namespace integration;

TEST_CASE("FRT oracle builds a valid flat HST",
          "[integration][tree][oracle][frt][flat-hst]")
{
    auto graph = makeCycleGraph();
    auto distances = unitDistances(*graph);

    FRT<FlatHST> oracle(*graph);
    auto tree = oracle.getTree(distances);

    requireValidFlatHST(tree, graph->getNumNodes());
}

TEST_CASE("FastCKR oracle builds a valid flat HST",
          "[integration][tree][oracle][ckr][flat-hst]")
{
    auto graph = makeCycleGraph();
    auto distances = unitDistances(*graph);

    FastCKR<FlatHST> oracle(*graph);
    auto tree = oracle.getTree(distances);

    requireValidFlatHST(tree, graph->getNumNodes());
}

TEST_CASE("MST oracle builds a valid flat HST",
          "[integration][tree][oracle][mst][flat-hst]")
{
    auto graph = makeCycleGraph();
    auto distances = unitDistances(*graph);

    TreeMST<FlatHST> oracle(*graph);
    auto tree = oracle.getTree(distances);

    requireValidFlatHST(tree, graph->getNumNodes());
}

TEST_CASE("FRT oracle builds a valid pointer HST",
          "[integration][tree][oracle][frt][pointer-hst]")
{
    auto graph = makeCycleGraph();
    auto distances = unitDistances(*graph);

    FRT<std::shared_ptr<HSTNode>> oracle(*graph);
    auto tree = oracle.getTree(distances);

    requireValidPointerHST(tree, graph->getNumNodes());
}

TEST_CASE("FastCKR oracle builds a valid pointer HST",
          "[integration][tree][oracle][ckr][pointer-hst]")
{
    auto graph = makeCycleGraph();
    auto distances = unitDistances(*graph);

    FastCKR<std::shared_ptr<HSTNode>> oracle(*graph);
    auto tree = oracle.getTree(distances);

    requireValidPointerHST(tree, graph->getNumNodes());
}

TEST_CASE("MST oracle builds a valid pointer HST",
          "[integration][tree][oracle][mst][pointer-hst]")
{
    auto graph = makeCycleGraph();
    auto distances = unitDistances(*graph);

    TreeMST<std::shared_ptr<HSTNode>> oracle(*graph);
    auto tree = oracle.getTree(distances);

    requireValidPointerHST(tree, graph->getNumNodes());
}