//
// Created by Mert Biyikli on 09.06.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_ENGINE_H
#define OBLIVIOUSROUTING_ROUTING_ENGINE_H

#include "routing_result.h"
#include "../core/utils.h"
#include <optional>
#include "../algorithms/oblivious/mwu/electrical_mwu.h"
#include "../algorithms/oblivious/mwu/tree_mwu.h"
#include "io/parse_argurment_io.h"


class RoutingEngine {
public:
    std::optional<IRoutingResult> solve(
        IGraph& graph,
        const Config& cfg,
        SolverType type
    );
};






#endif //OBLIVIOUSROUTING_ROUTING_ENGINE_H