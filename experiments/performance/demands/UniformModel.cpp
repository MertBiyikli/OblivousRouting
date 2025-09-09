//
// Created by Mert Biyikli on 08.09.25.
//

#include "UniformModel.h"
#include <random>

DemandMap UniformModel::generate(Graph& g, std::vector<std::pair<int, int>>& demands, double margin) {
    DemandMap demand2flow;

    std::uint64_t seed = std::random_device{}();
    std::mt19937_64 rng(seed);

    std::uniform_int_distribution<int> uniform_int(0, 401);
    std::uniform_real_distribution<double> uniform_real(0, 1);
    for(const auto& d : demands) {
        double flow = uniform_int(rng)+100;
        if (flow > 0) {
            const double flow_min = flow / margin;
            const double flow_max = flow * margin;
            const double r = uniform_real(rng);
            flow = r * (flow_max - flow_min) + flow_min;

        }

        demand2flow[d]=flow;
    }
    return demand2flow;
}