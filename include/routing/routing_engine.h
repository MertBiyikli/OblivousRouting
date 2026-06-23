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

class RoutingEngine
{
public:
    std::optional<RoutingRunResult> solve(IGraph& graph,const Config& cfg,const SolverType& type);
private:

    static bool isSemiObliviousSolver(SolverType type);
    std::optional<RoutingRunResult> solveSemiOblivious(IGraph& graph,const Config& cfg,SolverType type);


    static std::shared_ptr<IRoutingEngine> makeSemiRoutingEngine(SolverType type,IGraph& graph) {
        switch (type) {
            case SolverType::SEMI_ELECTRICAL:
                return std::make_shared<ExistingSolverRoutingEngine>(std::make_shared<ElectricalMWU>(graph, 0, true));

            case SolverType::SEMI_TREE:
                return std::make_shared<ExistingSolverRoutingEngine>(std::make_shared<TreeMWU<FlatHST>>(graph,0, std::make_unique<FastCKR<FlatHST>>(graph)));

            default:
                throw std::invalid_argument(
                    "Requested semi-oblivious routing engine for non-semi solver"
                );
        }
    }
};



#endif //OBLIVIOUSROUTING_ROUTING_ENGINE_H