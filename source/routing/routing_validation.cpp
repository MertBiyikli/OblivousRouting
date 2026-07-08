
#include "routing/routing_validation.h"

Result<void> RoutingValidation::input(const IGraph& _graph, const Config& _config) const {
    auto validate_graph = validateGraph(_graph);
    if (!validate_graph) {
        return getError(validate_graph);
    }

    auto validate_config = validateConfig(_config);
    if (!validate_config) {
        return getError(validate_config);
    }

    if (_config.evaluate_demand_models) {
        for (auto& [nam, demand]: _config.demand_maps) {
            auto validate_demands = validateDemand(demand, _graph);
            if (!validate_demands) {
                return getError(validate_demands);
            }
        }
    }


    return{};
}



Result<void> RoutingValidation::output(const IRoutingResult& _routing_results) const {
    auto validate_status = validateStatus(_routing_results);
    if (!validate_status) {
        return getError(validate_status);
    }

    auto validate_time_stats = validateTimeStats(_routing_results);
    if (!validate_time_stats) {
        return getError(validate_time_stats);
    }

    auto validate_routing = validateRouting(_routing_results);
    if (!validate_routing) {
        return getError(validate_routing);
    }

    return{};
}


Result<void> RoutingValidation::validateGraph(const IGraph &_graph) const {
    if (_graph.getNumNodes() == 0) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "Graph has no nodes.");
    }

    if (_graph.getNumDirectedEdges() == 0 || _graph.getNumUndirectedEdges() == 0) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "Graph has no edges.");
    }

    for (int e = 0; e < _graph.getNumDirectedEdges(); ++e) {
        const auto [head, tail] = _graph.getEdgeEndpoints(e);

        if (head < 0 || head >= _graph.getNumNodes()) {
            return makeErrorMessage(ErrorCode::InvalidGraph, "Edge has invalid head node.");
        }

        if (tail < 0 || tail >= _graph.getNumNodes()) {
            return makeErrorMessage(ErrorCode::InvalidGraph, "Edge has invalid tail node.");
        }

        if (!std::isfinite(_graph.getEdgeCapacity(e)) || _graph.getEdgeCapacity(e) <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidGraph, "Edge has non-positive or non-finite capacity.");
        }

        if (head == tail) {
            return makeErrorMessage(ErrorCode::InvalidGraph, "Graph contains a self-loop.");
        }
    }

    return{};
}

Result<void> RoutingValidation::validateDemand(const demands &_demands, const IGraph& graph) const {
    if (_demands.size()==0) {
        return std::unexpected(Error{ErrorCode::InvalidDemand, "Demand set is empty."});
    }

    for (int i = 0; i<_demands.size(); i++) {
        const auto& [source, target] = _demands.getDemandPair(i);
        double value = _demands.getDemandValue(i);

        if (source < 0 || source >= graph.getNumNodes()) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Demand has invalid source node.");
        }

        if (target < 0 || target >= graph.getNumNodes()) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Demand has invalid target node.");
        }

        if (source == target) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Demand source equals target.");
        }

        if (!std::isfinite(value) || value < 0.0) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Demand has negative or non-finite value.");
        }
    }

    return{};
}

Result<void> RoutingValidation::validateConfig(const Config &_config) const {

    if (_config.num_threads <= 0) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Number of threads must be positive.");
    }

    if (_config.evaluate_demand_models) {
        if (_config.demand_models.empty()) {
            return makeErrorMessage(ErrorCode::InvalidDemand,"Demand models are empty, but solver expects to route a demand.");
        }

        if (_config.offline_opt_per_model.empty()) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Offline optima are empty, but solver expects to route a demand.");
        }
    }


    for ( auto& [model_name, opt_value] : _config.offline_opt_per_model) {
        if (model_name.empty()) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Offline optimum contains an empty model name.");
        }

        if (!std::isfinite(opt_value) || opt_value < 0.0) {
            return makeErrorMessage(ErrorCode::RuntimeError,"Offline optimum for model '" + model_name + "' is negative or non-finite.");
        }

        if (!_config.demand_maps.contains(model_name)) {
            return makeErrorMessage(ErrorCode::RuntimeError,"Offline optimum was provided for model without matching demand map: " + model_name);
        }
    }

    for (const auto& [model_name, demand_map] : _config.demand_maps) {
        if (model_name.empty()) {
            return makeErrorMessage(ErrorCode::InvalidDemand, "Demand map contains an empty model name.");
        }

        if (demand_map.size() == 0) {
            return makeErrorMessage(ErrorCode::InvalidDemand,"Demand map for model '" + model_name + "' is empty.");
        }
    }
    return{};
}

Result<void> RoutingValidation::validateStatus(const IRoutingResult &result) const {
    if (result.status != ResultStatus::OK) {
        return makeErrorMessage(ErrorCode::SolverFailed, "Solver status is non-OK.");
    }

    if (result.oblivious_ratio >= 0.0 && !std::isfinite(result.oblivious_ratio)) {
        return makeErrorMessage(ErrorCode::SolverFailed, "Result has non-finite oblivious ratio.");
    }

    return {};
}

Result<void> RoutingValidation::validateTimeStats(const IRoutingResult& result) const {
    if (result.total_runtime_microseconds < 0 ||
        result.solve_runtime_microseconds < 0 ) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Runtime statistics contain negative values.");
        }

    if (!result.mwu_metrics.empty()) {
        if (result.mwu_metrics.getIterationCount() < 0) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Iteration count is negative.");
        }

        if (result.mwu_metrics.averageOracleTime() == -1.0) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Negative oracle running time.");
        }
    }
    return {};
}

Result<void> RoutingValidation::validateRouting(const IRoutingResult &result) const {
    if (result.preprocessing_runtime_microseconds > 0.0) {
        // Semi Oblivious routing:

    }else {
        if (!(result.scheme
            && result.scheme->isValid())){
            return makeErrorMessage(ErrorCode::InvalidRouting, "Routing table is invalid.");
            }
    }
    return {};
}