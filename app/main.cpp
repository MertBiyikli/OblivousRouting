#include "../include/routing/routing_engine.h"

int main(int argc, char **argv) {

    // Parse command line arguments
    RoutingEngine engine;
    auto res = engine.entry(argc, argv);

    if (!res) {
        std::cerr << "Error: " << getError(res).error().message << std::endl;
        return 1;
    }
    return 0;
}
