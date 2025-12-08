//
// Created by Mert Biyikli on 12.08.25.
//

#include "GravityModel.h"
#include "../../src/datastructures/IGraph.h"
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <random>

DemandMap GravityModel::generate(IGraph& g, std::vector<std::pair<int, int>>& demands, double margin) {
    DemandMap demand2flow;

    std::unordered_set<int> nodes(g.getNumNodes());
    for(const auto& d : demands) {
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
        }else{
            if(debug)
                std::cout << "Node " << i << " has capacity: " << nodeToCapacity[i] << "\n";
        }
    }

    std::uint64_t seed = std::random_device{}();
    std::mt19937_64 rng(seed);

    double capU, capV;
    double lo, hi;
    for(const auto& d : demands) {
        capU = nodeToCapacity[d.first];
        capV = nodeToCapacity[d.second];

        // Match your Java code: product / sumOfCapacities
        double baseFlow = (capU * capV) / sumCapacity;

        lo = baseFlow / margin;
        hi = baseFlow * margin;

        std::uniform_real_distribution<double> dist(0.0, 1.0);
        if(debug) {
            std::cout << "Generating demand for edge (" << d.first << ", " << d.second << ") with base flow: "
                      << baseFlow << ", lo: " << lo << ", hi: " << hi << "\n";
        }
        const double flow = dist(rng)*(hi-lo)+lo;
        // TODO: change this back

        demand2flow.addDemand(d.first, d.second, flow);
        //demand2flow[{d.first, d.second}] = flow; // Add the demand to the edge
    }
    return demand2flow;
}