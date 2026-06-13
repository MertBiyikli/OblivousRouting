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

struct SemiObliviousRoutingResult {
    std::unique_ptr<RoutingScheme> scheme;
    DemandModelType demand_type{};
    double congestion = -1.0;
    double runtime_microseconds = -1.0;

    std::size_t candidate_paths = 0;
    double average_paths_per_pair = 0.0;
    std::string path_selection_strategy;
};



class SemiObliviousRoutingSolver {
public:
    SemiObliviousRoutingSolver(
        std::shared_ptr<IRoutingEngine> routingEngine,
        std::shared_ptr<ISemiObliviousRoutingLoadOptimizer> loadOptimizer
    )
        : routingEngine_(std::move(routingEngine)),
          loadOptimizer_(std::move(loadOptimizer)) {
        if (!routingEngine_ || !loadOptimizer_) {
            throw std::invalid_argument("SemiObliviousRoutingSolver received null dependency");
        }
    }

    CandidateRoutingScheme preprocess(const IGraph& graph) {
        graph_ = &graph;
        candidateScheme_ = routingEngine_->preprocess(graph);
        return candidateScheme_.value();
    }

    SemiObliviousRoutingResult route(const demands& demand, DemandModelType demandType) const {
        if (!graph_ || !candidateScheme_) {
            throw std::logic_error("SemiObliviousRoutingSolver::preprocess must be called before route");
        }

        SemiObliviousRoutingResult res;
        res.demand_type = demandType;
        res.path_selection_strategy = routingEngine_->getSolverBase();

        res.candidate_paths = candidateScheme_->numPaths();
        res.average_paths_per_pair = candidateScheme_->averagePathsPerPair(
            graph_->getNumNodes()
        );

        const auto start = timeNow();

        auto optResult = loadOptimizer_->optimize(
            *graph_,
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
    const IGraph* graph_ = nullptr;

    std::shared_ptr<IRoutingEngine> routingEngine_;
    std::shared_ptr<ISemiObliviousRoutingLoadOptimizer> loadOptimizer_;

    std::optional<CandidateRoutingScheme> candidateScheme_;
};

inline void printSemiObliviousResult(
    const SemiObliviousRoutingResult& r
) {
    std::cout << "Routing base: " << r.path_selection_strategy << std::endl;
    std::cout << "Demand [" << demandModelName(r.demand_type) << "]\n";
    std::cout << "  Congestion: " << r.congestion << '\n';
    std::cout << "  Runtime: " << r.runtime_microseconds << " us\n";
    std::cout << "  Candidate paths: " << r.candidate_paths << '\n';
    std::cout << "  Avg paths/pair: " << r.average_paths_per_pair << '\n';
}

#endif //OBLIVIOUSROUTING_SEMI_OBLIVIOUS_SOLVER_H