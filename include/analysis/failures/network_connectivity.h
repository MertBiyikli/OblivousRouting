//
// Created by Mert on 28.10.26.
//

#ifndef E_ROUTING_NETWORK_CONNECTIVITY_H
#define E_ROUTING_NETWORK_CONNECTIVITY_H

#include "data_structures/graph/graph.h"

#include <queue>
#include <vector>

template <typename GraphType>
std::vector<int> connectedComponents(const GraphType& graph) {
    const int n = graph.getNumNodes();

    std::vector<std::vector<int>> adjacency(n);

    // Stable directed IDs are not necessarily contiguous
    // among active edges.
    const int bound = graph.globalEdgeCount() * 2;

    for (int e = 0; e < bound; e += 2) {
        if (!graph.edgeAlive(e))
            continue;

        const auto [u, v] = graph.getEdgeEndpoints(e);

        adjacency[u].push_back(v);
        adjacency[v].push_back(u);
    }

    std::vector<int> component(n, -1);
    std::queue<int> queue;

    int nextComponent = 0;

    for (int start = 0; start < n; ++start) {
        if (component[start] != -1)
            continue;

        component[start] = nextComponent;
        queue.push(start);

        while (!queue.empty()) {
            const int u = queue.front();
            queue.pop();

            for (const int v : adjacency[u]) {
                if (component[v] != -1)
                    continue;

                component[v] = nextComponent;
                queue.push(v);
            }
        }

        ++nextComponent;
    }

    return component;
}

inline bool isConnected(const optimized::Graph<EdgeData>& graph) {
    const auto components = connectedComponents(graph);

    for (const int component : components) {
        if (component != 0)
            return false;
    }

    return true;
}

inline Result<void> verifyDemandConnectivity(const optimized::Graph<EdgeData>& graph, const demands& demand_map) {
    const auto component = connectedComponents(graph);

    const auto& sources = demand_map.source;
    const auto& targets = demand_map.target;
    const auto& values = demand_map.demand_values;

    if (sources.size() != targets.size() || sources.size() != values.size()) {
        return makeErrorMessage(ErrorCode::InvalidRouting, "Malformed demand vectors.");
    }

    for (std::size_t i = 0; i < values.size(); ++i) {
        const int s = sources[i];
        const int t = targets[i];
        const double volume = values[i];

        if (s < 0 || t < 0 || s >= graph.getNumNodes() || t >= graph.getNumNodes()) {
            return makeErrorMessage(ErrorCode::InvalidRouting, "Demand endpoint outside graph.");
        }

        if (!std::isfinite(volume) || volume < 0.0) {
            return makeErrorMessage(ErrorCode::InvalidRouting, "Invalid demand volume.");
        }

        if (volume > 0.0 && component[s] != component[t]) {
            return makeErrorMessage(ErrorCode::InvalidRouting, "Link failure disconnects positive demand "
                                                               "between nodes " +
                                                                   std::to_string(s) + " and " + std::to_string(t) + ".");
        }
    }

    return {};
}

#endif // E_ROUTING_NETWORK_CONNECTIVITY_H
