//
// Created by Mert Biyikli on 11.07.26.
//
#include "algorithms/semi_oblivious/expander_hierarchy/preprocessing/hierarchy_preprocessor.h"

#include "xcut/core/config.hpp"
#include "xcut/data_structures/graph.hpp"
#include "xcut/expanders/expander_hierarchy.hpp"
#include "xcut/expanders/sparsifier.hpp"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <algorithm>
#include <limits>
#include <numeric>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace {

    Result<void>
    validateHierarchy(const optimized::Graph<EdgeData>& graph,const HierarchyResult& hierarchy) {
        if (hierarchy.levels.empty()) {
            return makeErrorMessage(ErrorCode::InvalidArgument,"Hierarchy is empty.");
        }

        for (const auto& level : hierarchy.levels) {
            std::vector<int> occurrence(graph.getNumNodes(),0);

            for (const auto& cluster : level.clusters) {
                if (cluster.original_vertices.empty()) {
                    return makeErrorMessage(ErrorCode::InvalidArgument,"Hierarchy contains an empty cluster.");
                }

                for (const int vertex : cluster.original_vertices) {
                    if (vertex < 0 ||
                        vertex >= graph.getNumNodes()) {
                            return makeErrorMessage(ErrorCode::InvalidGraph,"Hierarchy contains invalid vertex.");
                        }
                    ++occurrence[vertex];
                     }
            }

            for (int vertex = 0;vertex < graph.getNumNodes();++vertex) {
                if (occurrence[vertex] != 1) {
                    return makeErrorMessage(ErrorCode::InvalidGraph,
                        "Vertex " +
                            std::to_string(vertex) +
                            " occurs " +
                            std::to_string(occurrence[vertex]) +
                            " times at level " +
                            std::to_string(level.level)
                    );
                }
                 }
        }

        return {};
    }

        bool samePartition(const HierarchyLevel& lhs,const HierarchyLevel& rhs) {
            if (lhs.clusters.size() != rhs.clusters.size()) {
                return false;
            }

            auto canonicalize = [](const HierarchyLevel& level) {
                std::vector<std::vector<int>> parts;
                parts.reserve(level.clusters.size());

                for (const auto& cluster : level.clusters) {
                    auto vertices = cluster.original_vertices;
                    std::sort(vertices.begin(), vertices.end());
                    parts.push_back(std::move(vertices));
                }

                std::sort(parts.begin(), parts.end());
                return parts;
            };

            return canonicalize(lhs) == canonicalize(rhs);
        }

        void removeDuplicateTerminalLevels(HierarchyResult& hierarchy) {
            while (hierarchy.levels.size() >= 2) {
                const auto& previous = hierarchy.levels[hierarchy.levels.size() - 2];

                const auto& last = hierarchy.levels.back();

                if (!samePartition(previous, last)) {
                    break;
                }

                hierarchy.levels.pop_back();
            }

            for (int level_index = 0; level_index <static_cast<int>(hierarchy.levels.size()); ++level_index) {
                hierarchy.levels[level_index].level = level_index;

                for (auto& cluster : hierarchy.levels[level_index].clusters) {
                    cluster.level = level_index;
                }
            }
        }


    int ancestorAtLevel(Sparsifier& sparsifier,int original_vertex,int target_level) {
        int current = original_vertex;

        for (int level = 0; level < target_level; ++level) {
            current = sparsifier.parent(level, current);
        }

        return current;
    }

        void ensureXCutLogger() {
            if (!spdlog::get("xcut")) {
                auto logger = spdlog::stdout_color_mt("xcut");

                // Use info while debugging the hierarchy construction.
                logger->set_level(spdlog::level::off);

                // Optional:
                logger->set_pattern("[%n] [%l] %v");
            }
        }



bool containsAll(const std::vector<int>& parent,const std::vector<int>& children) {
    std::unordered_set<int> parent_set(parent.begin(), parent.end());

    return std::all_of(children.begin(),children.end(),[&](const int vertex) {
            return parent_set.contains(vertex);
        }
    );
}

} // namespace

std::pair<std::vector<std::pair<unsigned int,unsigned int>>, std::vector<double>> XCutHierarchyPreprocessor::toXCutEdges(const optimized::Graph<EdgeData>& graph) {
    std::vector<std::pair<unsigned int,unsigned int>> edges;
    std::vector<double> weight;
    edges.reserve(graph.getNumUndirectedEdges());
    weight.reserve(graph.getNumUndirectedEdges());

    for (int edge_id = 0; edge_id < graph.getNumDirectedEdges(); ++edge_id) {
        auto [u, v] = graph.getEdgeEndpoints(edge_id);

        if (u >= v) {
            continue;
        }

        edges.emplace_back(u, v);
        weight.emplace_back(graph.edgeData(edge_id).weight);
    }
    return {edges, weight};
}

std::vector<int> XCutHierarchyPreprocessor::computeInducedEdges(const optimized::Graph<EdgeData>& graph,const std::vector<int>& vertices) {
    std::vector<char> inside(graph.getNumNodes(), false);

    for (const int vertex : vertices) {
        if (vertex < 0 || vertex >= graph.getNumNodes()) {
            return {};
        }

        inside[vertex] = true;
    }

    std::vector<int> edges;
    edges.reserve(graph.getNumUndirectedEdges());

    /*
     * Iterate over actual directed edge IDs and retain one canonical
     * orientation for every physical undirected edge.
     */
    for (int e = 0;e < graph.getNumDirectedEdges();++e) {
        const auto [u, v] = graph.getEdgeEndpoints(e);

        if (u >= v) {
            continue;
        }

        if (inside[u] && inside[v]) {
            edges.push_back(e);
        }
    }

    return edges;
}

void XCutHierarchyPreprocessor::normalizeLevelOrder(HierarchyResult& hierarchy) {

    std::stable_sort(
        hierarchy.levels.begin(),
        hierarchy.levels.end(),
        [](const HierarchyLevel& lhs, const HierarchyLevel& rhs) {
            return lhs.clusters.size() > rhs.clusters.size();
        }
    );

    for (int level_index = 0;level_index < static_cast<int>(hierarchy.levels.size());++level_index) {
        hierarchy.levels[level_index].level = level_index;

        for (auto& cluster : hierarchy.levels[level_index].clusters) {
            cluster.level = level_index;
        }
    }
}

Result<void> XCutHierarchyPreprocessor::buildParentChildRelations(
    HierarchyResult& hierarchy
) {
    if (hierarchy.levels.empty()) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "Hierarchy has no levels.");
    }

    for (std::size_t level_index = 0;level_index + 1 < hierarchy.levels.size();++level_index) {

        auto& child_level = hierarchy.levels[level_index];
        auto& parent_level = hierarchy.levels[level_index + 1];

        for (auto& child : child_level.clusters) {
            HierarchyCluster* best_parent = nullptr;
            std::size_t best_parent_size = std::numeric_limits<std::size_t>::max();

            for (auto& candidate : parent_level.clusters) {
                if (!containsAll(
                        candidate.original_vertices,
                        child.original_vertices)) {
                    continue;
                }

                if (candidate.original_vertices.size() < best_parent_size) {
                    best_parent = &candidate;
                    best_parent_size =
                        candidate.original_vertices.size();
                }
            }

            if (best_parent == nullptr) {
                return makeErrorMessage(ErrorCode::InvalidGraph, "No parent found for child cluster " + std::to_string(child.id));
            }

            child.parent = best_parent->id;
            best_parent->children.push_back(child.id);
        }
    }

    auto& root_level = hierarchy.levels.back();

    if (root_level.clusters.size() != 1) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "Root level must contain exactly one cluster, but found " + std::to_string(root_level.clusters.size()));
    }

    hierarchy.root = root_level.clusters.front().id;
    return {};
}

void XCutHierarchyPreprocessor::buildLookupStructures(const optimized::Graph<EdgeData>& graph,HierarchyResult& hierarchy) {
    hierarchy.cluster_location.clear();
    hierarchy.vertex_to_leaf.assign(graph.getNumNodes(), -1);

    for (int level_index = 0; level_index < static_cast<int>(hierarchy.levels.size()); ++level_index) {
        auto& level = hierarchy.levels[level_index];

        for (int cluster_index = 0; cluster_index < static_cast<int>(level.clusters.size()); ++cluster_index) {
            const auto& cluster = level.clusters[cluster_index];

            hierarchy.cluster_location.emplace(
                cluster.id,
                std::pair{level_index, cluster_index}
            );
        }
    }

    if (hierarchy.levels.empty()) {
        return;
    }

    for (const auto& leaf : hierarchy.levels.front().clusters) {
        for (const int vertex : leaf.original_vertices) {
            hierarchy.vertex_to_leaf[vertex] = leaf.id;
        }
    }
}

void XCutHierarchyPreprocessor::choosePortals(const optimized::Graph<EdgeData>& graph,HierarchyResult& hierarchy) {
    for (auto& level : hierarchy.levels) {
        for (auto& cluster : level.clusters) {
            int best_vertex = -1;
            double best_internal_capacity = -1.0;

            for (const int vertex : cluster.original_vertices) {
                double internal_capacity = 0.0;

                for (const int edge_id : cluster.induced_edges) {
                    const auto [u, v] = graph.getEdgeEndpoints(edge_id);

                    if (u == vertex || v == vertex) {
                        internal_capacity += graph.edgeData(edge_id).capacity;
                    }
                }

                if (internal_capacity > best_internal_capacity) {
                    best_internal_capacity = internal_capacity;
                    best_vertex = vertex;
                }
            }

            cluster.portal = best_vertex;
        }
    }
}

Result<HierarchyResult>
XCutHierarchyPreprocessor::build(const optimized::Graph<EdgeData>& graph) const {
    if (graph.getNumNodes() == 0 || graph.getNumUndirectedEdges() == 0) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "Graph must have at least one vertex and one edge.");
    }

    ensureXCutLogger();

    const auto xcut_edges = toXCutEdges(graph);

    Graph xcut_graph(xcut_edges.first, xcut_edges.second,false);

    if (xcut_graph.has_degree_zero()) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "Graph has vertices with degree zero, which is not allowed for XCut.");
    }

    XCUT::Config config(0);
    config.m_verbose = false;
    config.m_debug = false;

    Sparsifier sparsifier =
        expander_hierarchy(&xcut_graph, &config);

   HierarchyResult hierarchy;
    hierarchy.levels.reserve(sparsifier.size());

    int next_cluster_id = 0;

    for (int xcut_level = 0;
         xcut_level < static_cast<int>(sparsifier.size());
         ++xcut_level) {

        HierarchyLevel level;
        level.level = xcut_level;

        std::unordered_map<int, std::vector<int>> cluster_vertices;

        for (int original_vertex = 0;
             original_vertex < graph.getNumNodes();
             ++original_vertex) {

            int level_vertex =
                static_cast<int>(original_vertex);

            /*
             * Lift the original vertex to graph(xcut_level).
             *
             * parent(level, u) is simultaneously:
             *   1. the part ID of u at this level, and
             *   2. the vertex ID representing that part in the next graph.
             */
            for (int level = 0; level < xcut_level; ++level) {
                level_vertex = sparsifier.parent(
                    static_cast<int>(level),
                    level_vertex
                );
            }

            const Graph* level_graph =
                sparsifier.graph(static_cast<int>(xcut_level));

            if (level_vertex >= level_graph->size()) {
                return makeErrorMessage(
                    ErrorCode::InvalidGraph,
                    "Invalid contracted vertex " +
                        std::to_string(level_vertex) +
                        " at XCut level " +
                        std::to_string(xcut_level)
                );
            }

            /*
             * This is the expander-decomposition cluster.
             * Do not use sparsifier.clustering().
             */
            const int partition_id = sparsifier.parent(xcut_level,level_vertex);

            cluster_vertices[partition_id].push_back(original_vertex);
        }

        level.clusters.reserve(cluster_vertices.size());

        for (auto& [partition_id, vertices] : cluster_vertices) {
            std::sort(vertices.begin(), vertices.end());

            HierarchyCluster cluster;
            cluster.id = next_cluster_id++;
            cluster.level = xcut_level;
            cluster.xcut_partition_id =
                static_cast<int>(partition_id);
            cluster.original_vertices = std::move(vertices);
            cluster.induced_edges = computeInducedEdges(graph,cluster.original_vertices);

            level.clusters.push_back(std::move(cluster));
        }

        hierarchy.levels.push_back(std::move(level));
    }
    auto validation = validateHierarchy(graph, hierarchy);
    if (!validation) {
        return getError(validation);
    }

    removeDuplicateTerminalLevels(hierarchy);

    buildLookupStructures(graph, hierarchy);

    auto relations = buildParentChildRelations(hierarchy);
    if (!relations) {
        return getError(relations);
    }

    buildLookupStructures(graph, hierarchy);
    //choosePortals(graph, hierarchy);

    for (const int leaf_id : hierarchy.vertex_to_leaf) {
        if (leaf_id < 0) {
            return makeErrorMessage(ErrorCode::InvalidGraph, "Some vertices are not assigned to any leaf cluster.");
        }
    }

    return hierarchy;
}