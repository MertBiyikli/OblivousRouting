//
// Created by Mert Biyikli on 08.06.26.
//

#ifndef OBLIVIOUSROUTING_TREE_ORACLE_TEST_HELPERS_H
#define OBLIVIOUSROUTING_TREE_ORACLE_TEST_HELPERS_H

#pragma once

#include <algorithm>
#include <set>
#include <string>
#include <vector>

#include "data_structures/graph/graph.h"
#include "algorithms/oblivious/mwu/oracle/tree/frt/frt.h"
#include "algorithms/oblivious/mwu/oracle/tree/fast_ckr/fast_ckr.h"
#include "algorithms/oblivious/mwu/oracle/tree/mst/mst_oracle.h"
#include "catch2/catch_test_macros.hpp"

// -----------------------------------------------------------------------------
// Shared graph builders
//
// TreeOracle implementations (FRT, FastCKR, TreeMST) operate on
// optimized::Graph<EdgeData>. Build test graphs using that type directly.
// -----------------------------------------------------------------------------

using TestGraph = optimized::Graph<EdgeData>;

inline TestGraph makeGraphFromEdges(int n, const std::vector<TestGraph::InputEdge>& edges) {
    return TestGraph(n, edges);
}

inline TestGraph makeSimplePathGraph5() {
    std::vector<TestGraph::InputEdge> edges = {
        {0, 1, EdgeData{1.0, 1.0}},
        {1, 2, EdgeData{1.0, 1.5}},
        {2, 3, EdgeData{1.0, 1.0}},
        {3, 4, EdgeData{1.0, 2.0}},
    };
    return makeGraphFromEdges(5, edges);
}

inline TestGraph makePathGraph(int n) {
    std::vector<TestGraph::InputEdge> edges;
    edges.reserve(n > 0 ? n - 1 : 0);

    for (int i = 0; i + 1 < n; ++i) {
        edges.push_back({i, i + 1, EdgeData{1.0, 1.0}});
    }

    return makeGraphFromEdges(n, edges);
}

inline TestGraph makeGridGraph(int size) {
    std::vector<TestGraph::InputEdge> edges;

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            const int node = i * size + j;

            if (j + 1 < size) {
                const int right = i * size + (j + 1);
                edges.push_back({node, right, EdgeData{1.0, 1.0}});
            }

            if (i + 1 < size) {
                const int down = (i + 1) * size + j;
                edges.push_back({node, down, EdgeData{1.0, 1.0}});
            }
        }
    }

    return makeGraphFromEdges(size * size, edges);
}

inline TestGraph makeFullyConnectedSmallGraph() {
    std::vector<TestGraph::InputEdge> edges;

    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            edges.push_back({i, j, EdgeData{1.0, 1.0}});
        }
    }

    return makeGraphFromEdges(4, edges);
}

inline TestGraph makeWeightedCycleGraph6() {
    std::vector<TestGraph::InputEdge> edges = {
        {0, 1, EdgeData{1.0, 1.0}},
        {1, 2, EdgeData{1.0, 2.0}},
        {2, 3, EdgeData{1.0, 1.0}},
        {3, 4, EdgeData{1.0, 2.0}},
        {4, 5, EdgeData{1.0, 1.0}},
        {5, 0, EdgeData{1.0, 3.0}},
        {0, 3, EdgeData{1.0, 4.0}},
    };
    return makeGraphFromEdges(6, edges);
}

inline std::vector<int> identityPermutation(int n) {
    std::vector<int> permutation;
    permutation.reserve(n);

    for (int v = 0; v < n; ++v) {
        permutation.push_back(v);
    }

    return permutation;
}

inline std::vector<double> currentGraphDistances(const TestGraph& graph) {
    std::vector<double> distances;
    distances.reserve(graph.getNumDirectedEdges());

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        distances.push_back(graph.edgeData(e).weight);
    }

    return distances;
}

inline std::vector<double> unitDistances(const TestGraph& graph) {
    return std::vector<double>(graph.getNumDirectedEdges(), 1.0);
}

inline std::vector<double> increasingDistances(const TestGraph& graph) {
    std::vector<double> distances(graph.getNumDirectedEdges());

    // Each undirected edge shares a single weight value across both of its
    // directed ids (directedId >> 1 == undirected edge id), so both
    // directions must be assigned the same distance for the value to stick.
    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        const int undirected = e >> 1;
        distances[e] = 1.0 + static_cast<double>(undirected % 5);
    }

    return distances;
}

// -----------------------------------------------------------------------------
// Oracle test case wrappers
// -----------------------------------------------------------------------------


struct FRTOracleCase {
    using Oracle = FRT<FlatHST>;

    static std::string name() {
        return "FRT";
    }

    static constexpr bool supports_level_partition = true;
    static constexpr bool supports_mendel_scaling = true;
};

struct FastCKROracleCase {
    using Oracle = FastCKR<FlatHST>;

    static std::string name() {
        return "FastCKR";
    }

    static constexpr bool supports_level_partition = true;
    static constexpr bool supports_mendel_scaling = true;
};

struct MSTOracleCase {
    using Oracle = TreeMST<FlatHST>;

    static std::string name() {
        return "TreeMST";
    }

    // TreeMST bypasses computeLevelPartition(...) and overrides getTree(...).
    static constexpr bool supports_level_partition = false;

    // TreeMST currently has only TreeMST(IGraph&), no TreeMST(IGraph&, bool).
    static constexpr bool supports_mendel_scaling = false;
};


// -----------------------------------------------------------------------------
// Shared validation helpers
// -----------------------------------------------------------------------------

inline bool isValidNodeId(int node, int n) {
    return node >= 0 && node < n;
}

inline bool containsValue(const std::vector<int>& values, int x) {
    return std::find(values.begin(), values.end(), x) != values.end();
}

inline bool hasUniqueValues(const std::vector<int>& values) {
    std::set<int> seen(values.begin(), values.end());
    return seen.size() == values.size();
}

inline void requireCentersAreValidNodes(const HSTLevel& level, int n) {
    for (int center : level.centers) {
        REQUIRE(isValidNodeId(center, n));
    }
}

inline void requireCentersAreContainedInPermutation(
    const HSTLevel& level,
    const std::vector<int>& permutation
) {
    for (int center : level.centers) {
        REQUIRE(containsValue(permutation, center));
    }
}

inline void requireOwnersAreValidOrUnassigned(const HSTLevel& level, int n) {
    for (int owner : level.owner) {
        REQUIRE(owner >= -1);
        REQUIRE(owner < n);
    }
}

inline void requireAllOwnersAreValidNodes(const HSTLevel& level, int n) {
    for (int owner : level.owner) {
        REQUIRE(isValidNodeId(owner, n));
    }
}

inline void requireAssignedOwnersAreCenters(const HSTLevel& level) {
    for (int owner : level.owner) {
        if (owner == -1) {
            continue;
        }

        REQUIRE(containsValue(level.centers, owner));
    }
}

inline void requirePartitionShapeIsValid(
    const HSTLevel& level,
    int number_of_nodes,
    const std::vector<int>& permutation
) {
    REQUIRE(level.owner.size() == static_cast<size_t>(number_of_nodes));

    requireCentersAreValidNodes(level, number_of_nodes);
    requireCentersAreContainedInPermutation(level, permutation);
    requireOwnersAreValidOrUnassigned(level, number_of_nodes);
}
#endif //OBLIVIOUSROUTING_TREE_ORACLE_TEST_HELPERS_H