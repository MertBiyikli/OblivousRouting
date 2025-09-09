//
// Created by Mert Biyikli on 08.09.25.
//

#include "BimodalModel.h"
#include <random>

DemandMap BimodalModel::generate(Graph& g, std::vector<std::pair<int, int>>& demands, double margin) {
    DemandMap demand2flow;

    std::uint64_t seed = std::random_device{}();
    std::mt19937_64 rng(seed);

    for(const auto& d : demands) {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        std::normal_distribution<double> gaussian(0.0, 1.0);

        for (const auto& d : demands) {
            double flow;
            if (uniform(rng) < 0.95) {
                // Low demand
                flow = gaussian(rng) * 20.+10; // mean 0, stddev 5
            } else {
                // High demand
                flow =  gaussian(rng) * 20.+400; // mean 0, stddev 5
            }
            if (flow < 0) flow = 0.0;

            demand2flow[d]=flow;
        }
    }
    return demand2flow;
}