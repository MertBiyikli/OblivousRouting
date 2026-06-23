//
// Created by Mert Biyikli on 12.06.26.
//

#ifndef OBLIVIOUSROUTING_SEMI_OBLIVIOUS_SOLVER_H
#define OBLIVIOUSROUTING_SEMI_OBLIVIOUS_SOLVER_H
#include <memory>
#include <optional>
#include <stdexcept>

#include "semi_routing_engine.h"
#include "load_optimizer.h"
#include "../../io/demand_io.h"
#include "core/solver.h"
#include "semi_oblivious_result.h"



class SemiObliviousRoutingSolver : public ISolver {
public:
    SemiObliviousRoutingSolver(IGraph& graph, std::shared_ptr<IRoutingEngine> routingEngine, std::shared_ptr<ISemiObliviousRoutingLoadOptimizer> loadOptimizer)
        : ISolver(graph),
          routingEngine_(std::move(routingEngine)),
          loadOptimizer_(std::move(loadOptimizer)) {
        if (!routingEngine_ || !loadOptimizer_) {
            throw std::invalid_argument(
                "SemiObliviousRoutingSolver received null dependency"
            );
        }
    }

    void setDemand(const demands& demand, DemandModelType demandType) {
        demand_ = demand;
        demandType_ = demandType;
    }

    CandidateRoutingScheme preprocess() {
        candidateScheme_ = routingEngine_->preprocess(graph);
        return *candidateScheme_;
    }

    std::unique_ptr<RoutingScheme> solve() override {
        if (!demand_ || !demandType_) {
            throw std::logic_error(
                "SemiObliviousRoutingSolver::setDemand must be called before solve"
            );
        }

        if (!candidateScheme_) {
            candidateScheme_ = routingEngine_->preprocess(graph);
        }

        auto optResult = loadOptimizer_->optimize(
            graph,
            *candidateScheme_,
            *demand_
        );

        return std::move(optResult.scheme);
    }

    SemiObliviousRoutingResult route(const demands& demand,DemandModelType demandType) {
        setDemand(demand, demandType);

        if (!candidateScheme_) {
            candidateScheme_ = routingEngine_->preprocess(graph);
        }

        SemiObliviousRoutingResult res;
        res.demand_type = demandType;
        res.path_selection_strategy = routingEngine_->getSolverBase();

        res.candidate_paths = candidateScheme_->numPaths();
        res.average_paths_per_pair =
            candidateScheme_->averagePathsPerPair(graph.getNumNodes());

        const auto start = timeNow();

        auto optResult = loadOptimizer_->optimize(
            graph,
            *candidateScheme_,
            demand
        );

        res.runtime_microseconds = duration(timeNow() - start);
        res.scheme = std::move(optResult.scheme);

        std::vector<double> cong;
        res.scheme->routeDemands(cong, demand);

        res.congestion = 0.0;
        for (double c : cong) {
            res.congestion = std::max(res.congestion, c);
        }

        return res;
    }

private:
    std::shared_ptr<IRoutingEngine> routingEngine_;
    std::shared_ptr<ISemiObliviousRoutingLoadOptimizer> loadOptimizer_;

    std::optional<CandidateRoutingScheme> candidateScheme_;
    std::optional<demands> demand_;
    std::optional<DemandModelType> demandType_;
};



#endif //OBLIVIOUSROUTING_SEMI_OBLIVIOUS_SOLVER_H