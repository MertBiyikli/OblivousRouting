//
// Created by Mert Biyikli on 24.06.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_RUNNER_H
#define OBLIVIOUSROUTING_ROUTING_RUNNER_H

#include "analysis/scenario_engine.h"
#include "core/errors.h"
#include "io/parse_argument_io.h"
#include "routing_result.h"


class IRoutingExperimentRunner {
  public:
    virtual ~IRoutingExperimentRunner() = default;

    virtual Result<IRoutingResult> run(optimized::Graph<EdgeData>& graph, const Config& cfg, SolverType type) const = 0;
};

class DemandEvaluator {
  public:
    static Result<NetworkData> evaluate(const optimized::Graph<EdgeData>& graph, const IRoutingResult& routing, const demands& demand) {
        double total_volume = 0.0;

        for (std::size_t i = 0; i < demand.size(); ++i) {
            const double value = demand.getDemandValue(i);

            if (!std::isfinite(value) || value < 0.0) {
                return makeErrorMessage(ErrorCode::InvalidDemand, "Invalid demand value.");
            }

            total_volume += value;

            if (!std::isfinite(total_volume)) {
                return makeErrorMessage(ErrorCode::InvalidDemand, "Total demand overflow.");
            }
        }

        return evaluateRouting(graph, routing, demand, total_volume);
    }
};

#endif // OBLIVIOUSROUTING_ROUTING_RUNNER_H
