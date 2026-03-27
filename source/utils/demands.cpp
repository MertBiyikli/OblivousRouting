#include "../../include/utils/demands.h"
#include "../../include/data_structures/graph/Igraph.h"
#include <cassert>
#include <random>
#include <unordered_set>

void demands::addDemand(int s, int t, double demand) {
    assert(demand_values.size() == source.size() && demand_values.size() == target.size());
    source.push_back(s);
    target.push_back(t);
    demand_values.push_back(demand);
}


size_t demands::size() const {
    assert(demand_values.size() == source.size() && demand_values.size() == target.size());
    return demand_values.size();
}

std::pair<int, int> demands::getDemandPair(size_t idx) const {
    assert(idx <= source.size() && idx <= target.size());
    return {source[idx], target[idx]};
}

double demands::getDemandValue(size_t idx) const {
    assert(idx <= demand_values.size());
    return demand_values[idx];
}


demands BimodalModel::generate(IGraph& g, std::vector<std::pair<int, int>>& _demands, double margin) {
    demands demand2flow;

    std::uint64_t seed = std::random_device{}();
    std::mt19937_64 rng(seed);


    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double> gaussian(0.0, 1.0);

    for (const auto& d : _demands) {
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
    }

    return demand2flow;
}

demands UniformModel::generate(IGraph& g, std::vector<std::pair<int, int>>& _demands, double margin) {
    demands demand2flow;

    std::uint64_t seed = std::random_device{}();
    std::mt19937_64 rng(seed);

    std::uniform_int_distribution<int> uniform_int(0, 401);
    std::uniform_real_distribution<double> uniform_real(0, 1);
    for(const auto& d : _demands) {
        double flow = uniform_int(rng)+100;
        if (flow > 0) {
            const double flow_min = flow / margin;
            const double flow_max = flow * margin;
            const double r = uniform_real(rng);
            flow = r * (flow_max - flow_min) + flow_min;

        }

        demand2flow.addDemand(d.first, d.second, flow);
    }
    return demand2flow;
}


demands GravityModel::generate(IGraph& g, std::vector<std::pair<int, int>>& _demands, double margin) {
    demands demand2flow;

    std::unordered_set<int> nodes(g.getNumNodes());
    for(const auto& d : _demands) {
        nodes.insert(d.first);
        nodes.insert(d.second);
    }

    std::vector<double> nodeToCapacity(g.getNumNodes(), 0.0);
    double sumCapacity = 0.0;
    for(const auto& node : nodes) {
        double sum = 0.0;
        for(const auto& u : g.neighbors(node)) {
            sum += g.getEdgeCapacity(node, u);
        }
        sumCapacity += sum;
        nodeToCapacity[node] = sum;
    }

    for(int i = 0; i < g.getNumNodes(); ++i) {
        if(nodeToCapacity[i] == 0.0) {
            throw std::runtime_error("Node " + std::to_string(i) + " has zero capacity, cannot generate demands.");
        }
    }

    std::uint64_t seed = std::random_device{}();
    std::mt19937_64 rng(seed);

    double capU, capV;
    double lo, hi;
    for(const auto& d : _demands) {
        capU = nodeToCapacity[d.first];
        capV = nodeToCapacity[d.second];

        // Match your Java code: product / sumOfCapacities
        double baseFlow = (capU * capV) / sumCapacity;

        lo = baseFlow / margin;
        hi = baseFlow * margin;

        std::uniform_real_distribution<double> dist(0.0, 1.0);
        const double flow = dist(rng)*(hi-lo)+lo;


        demand2flow.addDemand(d.first, d.second, flow);
    }
    return demand2flow;
}

demands GaussianModel::generate(IGraph& g, std::vector<std::pair<int, int>>& _demands, double margin) {
    demands demand2flow;
    std::uint64_t seed = std::random_device{}();
    std::mt19937_64 rng(seed);

    std::uniform_int_distribution<int> uniform_int(0, 401);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::normal_distribution<double> gaussian(0.0, 401);

    for(const auto& d : _demands) {

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