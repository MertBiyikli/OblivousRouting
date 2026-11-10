//
// Created by Mert on 30.10.26.
//

#ifndef E_ROUTING_COMMAND_H
#define E_ROUTING_COMMAND_H

#include "analysis/analysis_request.h"
#include "analysis/scenario_results.h"
#include "core/config.h"

#include <vector>


enum class Action {
    Run,
    AnalyzeN1,
    AnalyzeTraffic,
    AnalyzeScenario,
    Help,
    ListSolvers
};


struct CommodityScale {
    int source;
    int target;
    double factor;
};


struct Command {
    Action action = Action::Run;
    Config config{};
    ExecutionMode execution_mode = ExecutionMode::FAST;
    std::vector<CommodityScale> commodity_scales;
    std::optional<ScenarioAnalysisRequest> scenario_request;
};




#endif // E_ROUTING_COMMAND_H
