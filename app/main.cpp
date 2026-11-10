#include "../include/routing/routing_engine.h"
#include "analysis/resilience/resilience_result_json.h"
#include "application/resilience_application.h"
#include "application/routing_application.h"
#include "application/scenario_application.h"
#include "application/traffic_application.h"
#include "io/cli/help.h"
#include "io/cli/parser.h"
#include "io/result_io.h"
#include "routing/storage/routing_table_importer.h"
#include <string>

int main(int argc, char** argv) {
    auto command = parse(argc, argv);

    if (!command) {
        std::cerr << "Error: " << getError(command).error().message << std::endl;

        return 1;
    }

    switch (command->action) {
    case Action::Help:
        std::cout << help(argv[0]);
        return 0;

    case Action::ListSolvers:
        std::cout << solverList();
        return 0;

    case Action::AnalyzeScenario:
        return ScenarioApplication::run(command.value());

    case Action::AnalyzeN1:
        return ResilienceApplication::run(command.value());

    case Action::AnalyzeTraffic:
        return TrafficApplication::run(command.value());

    case Action::Run:
        return RoutingApplication::run(command.value());
    }

    std::cerr << "Unsupported CLI action.\n";
    return 1;
}
