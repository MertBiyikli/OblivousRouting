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

#include "data_structures/graph/graph_csr.h"
#include "algorithms/mwu/oracle/tree/frt/frt.h"
#include "algorithms/mwu/oracle/tree/fast_ckr/fast_ckr.h"
#include "algorithms/mwu/oracle/tree/mst/mst_oracle.h"
#include "catch2/catch_test_macros.hpp"

// -----------------------------------------------------------------------------
// Shared graph builders
// -----------------------------------------------------------------------------

inline GraphCSR makeSimplePathGraph5() {
    GraphCSR graph(5);

    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(1, 2, 1.0, 1.5);
    graph.addEdge(2, 3, 1.0, 1.0);
    graph.addEdge(3, 4, 1.0, 2.0);

    graph.finalize();
    return graph;
}

inline GraphCSR makePathGraph(int n) {
    GraphCSR graph(n);

    for (int i = 0; i + 1 < n; ++i) {
        graph.addEdge(i, i + 1, 1.0, 1.0);
    }

    graph.finalize();
    return graph;
}

inline GraphCSR makeGridGraph(int size) {
    GraphCSR graph(size * size);

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            const int node = i * size + j;

            if (j + 1 < size) {
                const int right = i * size + (j + 1);
                graph.addEdge(node, right, 1.0, 1.0);
            }

            if (i + 1 < size) {
                const int down = (i + 1) * size + j;
                graph.addEdge(node, down, 1.0, 1.0);
            }
        }
    }

    graph.finalize();
    return graph;
}

inline GraphCSR makeFullyConnectedSmallGraph() {
    GraphCSR graph(4);

    for (int i = 0; i < 4; ++i) {
        for (int j = i + 1; j < 4; ++j) {
            graph.addEdge(i, j, 1.0, 1.0);
        }
    }

    graph.finalize();
    return graph;
}

inline GraphCSR makeWeightedCycleGraph6() {
    GraphCSR graph(6);

    graph.addEdge(0, 1, 1.0, 1.0);
    graph.addEdge(1, 2, 1.0, 2.0);
    graph.addEdge(2, 3, 1.0, 1.0);
    graph.addEdge(3, 4, 1.0, 2.0);
    graph.addEdge(4, 5, 1.0, 1.0);
    graph.addEdge(5, 0, 1.0, 3.0);
    graph.addEdge(0, 3, 1.0, 4.0);

    graph.finalize();
    return graph;
}

inline std::vector<int> identityPermutation(int n) {
    std::vector<int> permutation;
    permutation.reserve(n);

    for (int v = 0; v < n; ++v) {
        permutation.push_back(v);
    }

    return permutation;
}

inline std::vector<double> currentGraphDistances(const IGraph& graph) {
    std::vector<double> distances;
    distances.reserve(graph.getNumDirectedEdges());

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        distances.push_back(graph.getEdgeDistance(e));
    }

    return distances;
}

inline std::vector<double> unitDistances(const IGraph& graph) {
    return std::vector<double>(graph.getNumDirectedEdges(), 1.0);
}

inline std::vector<double> increasingDistances(const IGraph& graph) {
    std::vector<double> distances(graph.getNumDirectedEdges());

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        distances[e] = 1.0 + static_cast<double>(e % 5);
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