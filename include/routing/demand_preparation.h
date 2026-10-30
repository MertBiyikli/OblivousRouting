//
// Created by Mert on 29.10.26.
//

#ifndef E_ROUTING_DEMAND_PREPARATION_H
#define E_ROUTING_DEMAND_PREPARATION_H

#include "core/errors.h"
#include "io/demand_io.h"

#include <cmath>
#include <utility>

inline Result<demands> prepareDemands(optimized::Graph<EdgeData>& graph, DemandModelType demandType, double scale = 1.0) {
    if (!std::isfinite(scale) || scale <= 0.0) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Demand scale must be finite and positive.");
    }

    auto pairs = generateAllDemandPairs(graph);

    auto model = makeDemandModel(demandType);

    auto generated = model->generate(graph, pairs);

    if (!generated) {
        return getError(generated);
    }

    demands scaled;

    for (std::size_t i = 0; i < generated->size(); ++i) {
        const auto [source, target] = generated->getDemandPair(i);

        const double value = generated->getDemandValue(i) * scale;

        if (!std::isfinite(value) || value < 0.0) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Invalid scaled demand.");
        }

        scaled.addDemand(source, target, value);
    }

    return scaled;
}

#endif // E_ROUTING_DEMAND_PREPARATION_H
