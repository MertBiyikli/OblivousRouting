//
// Created by Mert Biyikli on 09.06.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_ENGINE_H
#define OBLIVIOUSROUTING_ROUTING_ENGINE_H

#include "routing_result.h"
#include "../core/errors.h"
#include <optional>
#include "../algorithms/oblivious/mwu/electrical_mwu.h"
#include "../algorithms/oblivious/mwu/tree_mwu.h"
#include "io/parse_argument_io.h"


class RoutingEngine {
public:

    Result<void> entry(int argc, char **argv);


    Result<IRoutingResult> solve(optimized::Graph<EdgeData>& graph,const Config& cfg,SolverType type);

    Config cfg;
    std::unique_ptr<optimized::Graph<EdgeData>> m_graph;
    Result<std::unique_ptr<optimized::Graph<EdgeData>>> getGraph() {
        if (!m_graph) {
            return makeErrorMessage(ErrorCode::RuntimeError, "RoutingEngine::getGraph: graph not set");
        }else {
            return std::move(m_graph);
        }
    }


    Config getConfig() {return cfg;}
};






#endif //OBLIVIOUSROUTING_ROUTING_ENGINE_H