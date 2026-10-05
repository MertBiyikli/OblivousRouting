#include "../include/routing/routing_engine.h"
#include "../include/io/cli.h"

int main(int argc, char **argv) {

    // Parse command line arguments
    auto cmd = parse(argc, argv);
    if (!cmd) {
        std::cerr << "Error: " << getError(cmd).error().message << std::endl;
        return 1;
    }

    switch (cmd->action) {

    case Action::Help:
        std::cout << help(argv[0]);
        return 0;

    case Action::ListSolvers:
        std::cout << solverList();
        return 0;

    case Action::Run:
        break;
    }



    // Running engine
    RoutingEngine engine;
    auto res = engine.entry(cmd.value().config);
    if (!res) {
        std::cerr << "Error: " << getError(res).error().message << std::endl;
        return 1;
    }
    return 0;//something new
}
