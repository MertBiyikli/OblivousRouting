//
// Created by Mert Biyikli on 29.09.25.
//

#include "linear_oblivious_routing_ratio.h"

void LinearObliviousRatio::init(GraphADJ &g, const std::unordered_map<std::pair<int, int>, std::unordered_map<std::pair<int, int>, double, PairHash> , PairHash > &routing) {
    this->routing = routing;
    this->g = g;
    // Collect unique commodity pairs
    for (const auto& v : g.getVertices()) {
        for (const auto& u : g.neighbors(v)) {
            if (v < u) {
                demand_vector[{v, u}] = g.getEdgeCapacity(v, u);
            }
        }
    }
}

double LinearObliviousRatio::solve() {
    std::unordered_map<std::pair<int, int>, double, PairHash> congestion_per_edge;

    for ( const auto& [edge, flowMap]:routing) {
        for (const auto& [com, flow_value]  : flowMap) {
            if ( flow_value < 1e-15 ) continue; // ignore zero flows
            auto undirected_edge = (edge.first < edge.second ? edge : std::make_pair(edge.second, edge.first));

            auto it = demand_vector.find( com );
            if ( it == demand_vector.end() ) continue;
            double demand = it->second;

            if ( flow_value == 0 ) continue;
            if (!congestion_per_edge.contains(undirected_edge) )
                congestion_per_edge[undirected_edge] = 0;

            congestion_per_edge[undirected_edge] += std::abs(demand*flow_value);
        }
    }
    max_congestion = -std::numeric_limits<double>::infinity();
    for ( const auto& [edge, total_flow]:congestion_per_edge) {
        if ( total_flow == 0 ) continue;
        double capacity = g.getEdgeCapacity(edge.first, edge.second);
        double cong = total_flow / capacity;
        max_congestion = std::max(max_congestion, cong);
    }

    return max_congestion;
}
