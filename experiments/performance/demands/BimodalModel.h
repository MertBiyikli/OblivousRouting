//
// Created by Mert Biyikli on 08.09.25.
//

#ifndef OBLIVIOUSROUTING_BINMOALMODEL_H
#define OBLIVIOUSROUTING_BINMOALMODEL_H

#include <random>
#include <vector>
#include "DemandModel.h"

class BimodalModel : public DemandModel {
public:
    bool debug = false;
    BimodalModel() = default;
    // Generate a gravity model demand matrix
    virtual DemandMap generate(IGraph& g, std::vector<std::pair<int, int>>& demands, double margin = 1.0) {
        DemandMap demand2flow;

        std::uint64_t seed = std::random_device{}();
        std::mt19937_64 rng(seed);


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

            demand2flow.addDemand(d.first, d.second, flow);
            // demand2flow[d]=flow;
        }

        return demand2flow;
    }
};

#endif //OBLIVIOUSROUTING_BINMOALMODEL_H