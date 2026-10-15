//
// Created by Mert on 15.10.26.
//

#ifndef E_ROUTING_ROUTING_SCHEME_VERIFIER_H
#define E_ROUTING_ROUTING_SCHEME_VERIFIER_H
#include "core/errors.h"
#include "routing/routing_table.h"

#include <algorithm>
#include <vector>

class RoutingSchemeVerifier {
  public:
    static Result<void> verify(const optimized::Graph<EdgeData>& graph, const RoutingScheme& scheme, double tolerance = 1e-7) {
        const int n = graph.getNumNodes();
        const int m = graph.getNumDirectedEdges();

        for (int source = 0; source < n; ++source) {
            for (int target = source + 1; target < n; ++target) {
                std::vector<double> balance(n, 0.0);
                double absoluteFlow = 0.0;

                for (int edge = 0; edge < m; ++edge) {
                    const int reverse = graph.reverse(edge).id;

                    if (edge > reverse)
                        continue;

                    const auto [u, v] = graph.getEdgeEndpoints(edge);
                    const double flow = scheme.getFlow(edge, source, target);

                    if (!std::isfinite(flow)) {
                        return makeErrorMessage(ErrorCode::InvalidRouting, "Imported routing scheme contains non-finite flow");
                    }

                    balance[u] += flow;
                    balance[v] -= flow;
                    absoluteFlow += std::abs(flow);
                }

                const double scale = std::max(1.0, absoluteFlow);
                const double allowedError = tolerance * scale;

                for (int vertex = 0; vertex < n; ++vertex) {
                    const double expected = vertex == source ? 1.0 : (vertex == target ? -1.0 : 0.0);

                    if (std::abs(balance[vertex] - expected) > allowedError) {
                        return makeErrorMessage(ErrorCode::InvalidRouting, "Imported scheme violates flow conservation for commodity " + std::to_string(source) + "->" + std::to_string(target) +
                                                                               " at vertex " + std::to_string(vertex));
                    }
                }
            }
        }

        return {};
    }
};

#endif // E_ROUTING_ROUTING_SCHEME_VERIFIER_H
