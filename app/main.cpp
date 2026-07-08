#include "../include/io/parse_argument_io.h"
#include "../include/routing/routing_engine.h"

int main(int argc, char **argv) {
    // Parse command line arguments
    RoutingEngine engine;
    auto result = engine.entry(argc, argv);
    if (!result) {
        std::cerr << "Error: " << result.error().message << std::endl;
        return 1;
    }
    return 0;
}
