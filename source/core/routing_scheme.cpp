//
// Created by Mert Biyikli on 20.03.26.
//

#include "../../include/core/routing_scheme.h"
#include "../../include/utils/my_math.h"

double RoutingScheme::getMaxCongestion(const std::vector<double>& congestion) const {
    double max_cong = 0.0;
    for (const auto& cong : congestion) {
        if (cong > max_cong) {
            max_cong = cong;
        }
    }
    return max_cong;
}





void LinearRoutingScheme::printRoutingTable() const {
    routing_table.printFlows(g);
}

double LinearRoutingScheme::computeObliviousRatio() {
    // for the liner routing scheme, we can compute the oblivious ratio, by pushing for each
    // edge the capacity of along the edge points
    const int m = g.getNumDirectedEdges();
    demands worst_case_demands;
    for (int e = 0; e < m; ++e) {
        auto [u, v] = g.getEdgeEndpoints(e);
        if (u > v) continue; // only consider one orientation for undirected edges
        double capacity = g.getEdgeCapacity(e);
        worst_case_demands.addDemand(u,v,capacity);
    }

    std::vector<double> congestion;
    routeDemands(congestion, worst_case_demands);
    double max_ratio = getMaxCongestion(congestion);
    return max_ratio;
}

void LinearRoutingScheme::initRoutingTable() {
    routing_table.init(g);
}



// Since we store the oblivious routing as a linear routing table w.r.t. root_x,
// and we encode the flow orientation  as (absolute) flows along the direction of the
// undirected edge which is represented as two directed edges (u,v) and (v,u),
// we need to consider both orientations when querying the flow for an undirected edge.
//
// Example: Given an edge e=(u,v), assume we have flow f_e_st going from v to u for demand s→t,
// which is the reverse orientation of the original edge (u,v)
//
// we encoded this as:
// In the representation as an undirected graph this would be negative flow along edge (u,v)
// but in our linear routing table representation w.r.t. root node x, we have:
//
// a positive flow f_e'_st along the anti-edge e'=(v,u) for demand s→t which is then represented by the linearity as
// - f_e_st = f_e'_st = f_e'_sx - f_e'_tx
//
//
// which means that if we want to get the flow for an edge e that is an undirected edge,
// and encode the flow orientation correctly, which is the by the signs of its value we need to consider the anti-edge as well.
double LinearRoutingScheme::getFlow(int e, int s, int t) const {
    int e_orig = e;
    int anti_e = g.getAntiEdge(e);

    double flow_sx = routing_table.getFlow(e_orig, s);
    double flow_tx = routing_table.getFlow(e_orig, t);

    double flow_sx_anti = routing_table.getFlow(anti_e, s);
    double flow_tx_anti = routing_table.getFlow(anti_e, t);

    double total_flow_sx = flow_sx - flow_sx_anti;
    double total_flow_tx = flow_tx - flow_tx_anti;

    double total_flow = total_flow_sx - total_flow_tx;
    return ( std::abs(total_flow) < EPS ? 0.0 : total_flow );
}

void LinearRoutingScheme::addFlow(int e, int s, int t, double flow_sx) {
    routing_table.addFlow(e, s, flow_sx);
    routing_table.addFlow(e, t, -flow_sx);
}

void LinearRoutingScheme::routeDemands(std::vector<double>& congestion,
                      const demands& demands) const {
    const int m = g.getNumDirectedEdges();

    // Initialize congestion for undirected edges only
    std::vector<double> directed_congestion(m, 0.0);

    // Process all directed edges, accumulating flow into undirected edge entries
    for (int e = 0; e < m; ++e) {
        double flow = 0.0;
        for (int d = 0; d < demands.size(); ++d) {
            const auto& [s, t] = demands.getDemandPair(d);
            const double coeff = demands.getDemandValue(d);
            if (coeff == 0.0) continue;

            flow += coeff * std::abs(getFlow(e, s, t));
        }
        // Add flow (normalize by capacity)
        directed_congestion[e] += flow / g.getEdgeCapacity(e);
    }

    // map the directed congestion back to the undirected edges
    congestion.assign(m , 0.0); // reset congestion to store undirected congestion
    for (int e = 0; e < m; ++e) {
        const auto& [u, v] = g.getEdgeEndpoints(e);
        int undirected_idx = (u < v ? e : g.getAntiEdge(e));
        congestion[undirected_idx] += directed_congestion[e];
    }
}

bool LinearRoutingScheme::isValid() {
    return routing_table.isValid(g);
}


void AllPairRoutingScheme::printRoutingTable() const {
    routing_table.printFlows(g);
}


double AllPairRoutingScheme::getFlow(int e, int s, int t) const{
    int e_orig = e;
    return routing_table.getFlow(e_orig, s, t);
}

void AllPairRoutingScheme::addFlow(int e, int s, int t, double flow_sx) {
    routing_table.addFlow(e, s,t, flow_sx);
}

void AllPairRoutingScheme::routeDemands(std::vector<double>& congestion, const demands& demands) const {
    const int m = g.getNumDirectedEdges();
    const int num_undirected = m / 2;  // number of undirected edges

    // Initialize congestion for undirected edges only
    congestion.assign(num_undirected, 0.0);

    // Process all directed edges, accumulating flow into undirected edge entries
    for (int e = 0; e < m; ++e) {
        double flow = 0.0;
        for (int d = 0; d < demands.size(); ++d) {
            const auto& [s, t] = demands.getDemandPair(d);
            const double coeff = demands.getDemandValue(d);
            if (coeff == 0.0) continue;

            flow += coeff * std::abs(routing_table.getFlow(e, s, t));
        }

        // Map directed edge to undirected edge index
        // Map directed edge to undirected edge index
        const auto& [u, v] = g.getEdgeEndpoints(e);
        int undirected_idx = (u < v ? e : g.getAntiEdge(e));

        // Add flow (normalize by capacity)
        congestion[undirected_idx] += flow / g.getEdgeCapacity(e);
    }
}


bool AllPairRoutingScheme::isValid() {
    return routing_table.isValid(g);
}
