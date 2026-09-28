//
// Created by Mert Biyikli on 12.06.26.
//

#ifndef OBLIVIOUSROUTING_SEMI_OBLIVIOUS_SOLVER_H
#define OBLIVIOUSROUTING_SEMI_OBLIVIOUS_SOLVER_H
#include <memory>
#include <optional>
#include <stdexcept>

#include "semi_routing_engine.h"
#include "postprocessing/load_optimizer.h"
#include "../../io/demand_io.h"
#include "core/solver.h"
#include "routing/routing_result.h"



class SemiObliviousRoutingSolver : public ISolver {
public:
    SemiObliviousRoutingSolver(optimized::Graph<EdgeData>& graph, std::shared_ptr<SemiSolverRoutingEngine> routingEngine, std::shared_ptr<ISemiObliviousRoutingLoadOptimizer> loadOptimizer)
        : ISolver(graph),
          routingEngine_(std::move(routingEngine)),
          loadOptimizer_(std::move(loadOptimizer)), current_result() {


        if (!routingEngine_ || !loadOptimizer_) {
            throw std::invalid_argument(
                "SemiObliviousRoutingSolver received null dependency"
            );
        }
    }

    Result<std::unique_ptr<RoutingScheme>> solve() override;

    void setDemand(const demands& demand,  const DemandModelType& demandType);
    Result<CandidateRoutingScheme> preprocess();
    Result<SemiObliviousRoutingResult> route(const demands& demand, const DemandModelType& demandType);

private:
    std::shared_ptr<SemiSolverRoutingEngine> routingEngine_;
    std::shared_ptr<ISemiObliviousRoutingLoadOptimizer> loadOptimizer_;

    std::optional<CandidateRoutingScheme> candidateScheme_;
    std::optional<demands> demand_;
    std::optional<DemandModelType> demandType_;

    SemiObliviousRoutingResult current_result;
};





#endif //OBLIVIOUSROUTING_SEMI_OBLIVIOUS_SOLVER_H