//
// Created by Mert on 26.10.26.
//

#ifndef E_ROUTING_SCENARIO_DEMANDS_H
#define E_ROUTING_SCENARIO_DEMANDS_H

#include <cmath>
#include <cstddef>

#include "core/errors.h"
#include "utils/demands.h"

// Include the header that declares your actual
// makeDemandModel() and generateAllDemandPairs() functions.
// Do not include that header from utils/demands.h.

template <typename DemandModelType, typename GraphType>
Result<demands> generateScenarioDemands(GraphType graph, DemandModelType demand_type, double demand_scale) {
    if (!std::isfinite(demand_scale) || demand_scale <= 0.0) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Demand scale must be finite and positive.");
    }

    auto pairs = generateAllDemandPairs(graph);
    auto model = makeDemandModel(demand_type);

    auto generated = model->generate(graph, pairs);

    if (!generated) {
        return getError(generated);
    }

    demands scaled;
    for (std::size_t i = 0; i < generated->size(); ++i) {

        const auto [source, target] = generated->getDemandPair(i);

        const double value = generated->getDemandValue(i) * demand_scale;

        if (!std::isfinite(value) || value < 0.0) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Invalid scaled demand.");
        }

        scaled.addDemand(source, target, value);
    }

    return scaled;
}

#endif // E_ROUTING_SCENARIO_DEMANDS_H
