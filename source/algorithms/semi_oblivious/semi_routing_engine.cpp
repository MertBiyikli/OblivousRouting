//
// Created by Mert Biyikli on 24.06.26.
//

#include "algorithms/semi_oblivious/semi_routing_engine.h"

#include "algorithms/oblivious/mwu/electrical_mwu.h"
#include "algorithms/oblivious/mwu/tree_mwu.h"
#include "algorithms/semi_oblivious/expander_hierarchy/tree_sparsifier_solver.h"
#include "utils/my_math.h"

Result<void> validateResidualOrientation(const IGraph& graph,const int source,const int target,const std::vector<double>& residual) {
    std::vector<double> divergence(
        graph.getNumNodes(),
        0.0
    );

    for (int edge = 0;
         edge < graph.getNumDirectedEdges();
         ++edge) {
        const double flow = residual[edge];

        if (flow <= 0.0) {
            continue;
        }

        const auto [u, v] =
            graph.getEdgeEndpoints(edge);

        divergence[u] += flow;
        divergence[v] -= flow;
         }

    constexpr double tolerance = 1e-7;

    /*
     * Unit s -> t flow:
     *
     *     div(s) = +1
     *     div(t) = -1
     */
    if (divergence[source] < 1.0 - tolerance ||
        divergence[target] > -1.0 + tolerance) {
        return makeErrorMessage(
            ErrorCode::InvalidRouting,
            "Flow orientation/conservation is invalid for pair " +
                std::to_string(source) +
                " -> " +
                std::to_string(target) +
                ". source divergence=" +
                std::to_string(divergence[source]) +
                ", target divergence=" +
                std::to_string(divergence[target])
        );
        }

    return {};
}

Result<CandidateRoutingScheme> SemiSolverRoutingEngine::preprocess(const IGraph& graph) {
    CandidateRoutingScheme candidateScheme;

    auto scheme = solver_->solve();
    if (!scheme || !scheme.value()) {
        return getError(scheme);
    }
    for (int s = 0; s < graph.getNumNodes(); s++) {
        for (int t = 0; t<graph.getNumNodes(); t++) {
            if (s == solver_->getRootNode() || t == solver_->getRootNode() || s == t) {
                continue;
            }

            auto extracted = extractPath(graph, *(scheme.value()), s, t, candidateScheme);
            if (!extracted) {
                return getError(extracted);
            }
        }
    }

    return candidateScheme;
}

Result<void> SemiSolverRoutingEngine::extractPath(const IGraph& graph,const RoutingScheme& scheme,const int source,const int target,CandidateRoutingScheme& output) const {
    const int directed_edges =
        graph.getNumDirectedEdges();

    std::vector<double> residual(
        directed_edges,
        0.0
    );

    constexpr double decomposition_epsilon = 1e-12;

    /*
     * Process each undirected edge exactly once.
     */
    for (int edge = 0;
         edge < directed_edges;
         ++edge) {
        const int anti =
            graph.getAntiEdge(edge);

        if (anti == INVALID_EDGE_ID) {
            return makeErrorMessage(
                ErrorCode::InvalidGraph,
                "Graph edge has no anti-edge."
            );
        }

        if (edge > anti) {
            continue;
        }

        const double flow =
            scheme.getFlow(
                edge,
                source,
                target
            );

        if (!std::isfinite(flow)) {
            return makeErrorMessage(
                ErrorCode::InvalidRouting,
                "Routing scheme returned a non-finite flow."
            );
        }

        if (flow > decomposition_epsilon) {
            residual[edge] = flow;
        } else if (flow < -decomposition_epsilon) {
            residual[anti] = -flow;
        }
    }

    auto orientation_check =
        validateResidualOrientation(
            graph,
            source,
            target,
            residual
        );

    if (!orientation_check) {
        return getError(orientation_check);
    }

    while (true) {
        std::vector<int> path_edges;

        std::vector<bool> visited(
            graph.getNumNodes(),
            false
        );

        const bool found =
            dfsDecompose(
                graph,
                source,
                target,
                residual,
                path_edges,
                visited
            );

        if (!found) {
            break;
        }

        double bottleneck =
            std::numeric_limits<double>::infinity();

        for (const int edge : path_edges) {
            bottleneck = std::min(
                bottleneck,
                residual[edge]
            );
        }

        if (!std::isfinite(bottleneck) ||
            bottleneck <= decomposition_epsilon) {
            break;
        }

        for (const int edge : path_edges) {
            residual[edge] -= bottleneck;

            if (std::abs(residual[edge]) <=
                decomposition_epsilon) {
                residual[edge] = 0.0;
            }
        }

        Path path;
        path.source = source;
        path.target = target;
        path.edges = std::move(path_edges);
        path.weight = bottleneck;

        output.addPath(
            source,
            target,
            std::move(path)
        );
    }

    return {};
}

const std::string SemiSolverRoutingEngine::getSolverBase() const {
    if (dynamic_cast<ElectricalMWU*>(
            solver_.get())) {
        return "Electrical";
            }

    if (dynamic_cast<TreeMWU<FlatHST>*>(
            solver_.get())) {
        return "Tree";
            }

    if (dynamic_cast<
            ElectrifiedExpanderHierarchySolver*
        >(solver_.get())) {
        return "Electrified Expander Hierarchy";
        }

}

bool SemiSolverRoutingEngine::dfsDecompose(const IGraph& g,int u,int t,std::vector<double>& residual,std::vector<int>& currentPath,std::vector<bool>& visited) const {
    if (u == t) {
        return true;
    }

    visited[u] = true;

    for (int neig: g.neighbors(u)) { // adapt name if needed
        int e = g.getEdgeId(u, neig);
        if (e == INVALID_EDGE_ID) continue;
        if (residual[e] <= EPS) continue;

        const auto [a, b] = g.getEdgeEndpoints(e);

        if (a != u) continue;

        const int v = b;

        if (visited[v]) continue;

        currentPath.push_back(e);

        if (dfsDecompose(g, v, t, residual, currentPath, visited)) {
            return true;
        }

        currentPath.pop_back();
    }

    return false;
}