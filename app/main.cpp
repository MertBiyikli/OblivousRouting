#include "../include/routing/routing_engine.h"

int main(int argc, char **argv) {

    // Parse command line arguments
    auto cfg = parse_parameter(argc, argv);
    if (!cfg) {
        std::cerr << "Error: " << getError(cfg).error().message << std::endl;
        return 1;
    }


    // Running engine
    RoutingEngine engine;
    auto res = engine.entry(cfg.value());
    if (!res) {
        std::cerr << "Error: " << getError(res).error().message << std::endl;
        return 1;
    }
    return 0;//something new
}
