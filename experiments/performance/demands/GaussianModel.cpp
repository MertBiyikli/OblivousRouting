//
// Created by Mert Biyikli on 08.09.25.
//

#include "GaussianModel.h"

#include <iostream>
#include <random>
#include "../../src/datastructures/IGraph.h"

DemandMap GaussianModel::generate(IGraph& g, std::vector<std::pair<int, int>>& demands, double margin) {
    DemandMap demand2flow;

    std::uint64_t seed = std::random_device{}();
    std::mt19937_64 rng(seed);

    std::uniform_int_distribution<int> uniform_int(0, 401);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double> gaussian(0.0, 401);

    std::cout << "Number of nodes:" << g.getNumNodes() << std::endl;

    for(const auto& d : demands) {

        double flow = (uniform_int(rng) + 100) + gaussian(rng); // mean 100, stddev 50
        if (flow > 0) {
            // flow = new Random().nextDouble()*(flow*margin - flow/margin) + flow/margin;
            double flow_min = flow / margin;
            double flow_max = flow * margin;
            double randFactor = uniform(rng);
            flow = randFactor * (flow_max - flow_min) + flow_min;

            // Insert two entries like in Java code
            demand2flow.addDemand(d.first, d.second, flow);
        }
    }
    return demand2flow;
}