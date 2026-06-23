//
// Created by Mert Biyikli on 23.06.26.
//
#include "routing/storage/linear_routing_table.h"
#include <cassert>
#include "utils/my_math.h"
#include <iostream>

void LinearRoutingTable::init(const IGraph& g) {
    const int numEdges = g.getNumDirectedEdges();
    n = g.getNumNodes();
    src_ids.assign(numEdges, {});
    src_flows.assign(numEdges, {});
}


void LinearRoutingTable::addFlow(int e, int s, double flow_sx) {
    assert(e >= 0 && e < src_ids.size() && s >= 0 && s < n);
    auto& ids  = src_ids[e];
    auto& vals = src_flows[e];

    // keep the ids sorted by s for binary search later
    size_t len = ids.size();
    size_t lo = 0, hi = len;

    const size_t linear_bound = 8;
    for (int i = 0; i < std::min(len, linear_bound); ++i) {
        if (ids[i] == s) {
            vals[i] += flow_sx;
            return;
        }
    }

    while (lo < hi) {
        const size_t mid = (lo + hi) >> 1;
        const int mid_val = ids[mid];
        if (mid_val < s)
            lo = mid + 1;
        else
            hi = mid;
    }

    // if found, update
    if (lo < len && ids[lo] == s) {
        vals[lo] += flow_sx;
    }else {
        // insert new value at the new position
        ids.insert(ids.begin()+static_cast<long>(lo), s);
        vals.insert(vals.begin()+static_cast<long>(lo), flow_sx);
    }
}


const double LinearRoutingTable::getFlow(int e, int s) const {
    assert(e >= 0 && e < src_ids.size() && s >= 0 && s < n);
    const auto& ids  = src_ids[e];
    const auto& vals = src_flows[e];

    // small linear scan first (your trick)
    const int len = static_cast<int>(ids.size());
    const int linear_bound = 8;
    for (int i = 0; i < std::min(len, linear_bound); ++i) {
        if (ids[i] == s) {
            return vals[i];
        }
    }

    // if you keep ids sorted by s, you can do binary search here
    int lo = 0, hi = len;
    while (lo < hi) {
        const size_t mid = (lo + hi) >> 1;
        const int mid_val = ids[mid];
        if (mid_val < s)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo < len && ids[lo] == s) {
        return vals[lo];
    } else {
        return 0.0;
    }
}

const int LinearRoutingTable::getSize() const {
    assert(src_ids.size() == src_ids.size());
    return src_ids.size();
}

void LinearRoutingTable::eraseFlow(int e, int s) {
    // remove flow at edge e with commodity s
    assert(e >= 0 && e < src_ids.size() && s >= 0 && s < n);
    auto& ids  = src_ids[e];
    auto& vals = src_flows[e];

    int index_s = -1;

    // small linear scan first (your trick)
    const int len = static_cast<int>(ids.size());
    const int linear_bound = 8;
    for (int i = 0; i < std::min(len, linear_bound); ++i) {
        if (ids[i] == s) {
            index_s = i;
        }
    }

    if (index_s != -1) {
        ids.erase(ids.begin()+index_s);
        vals.erase(vals.begin()+index_s);
        return;
    }

    // if you keep ids sorted by s, you can do binary search here
    int lo = 0, hi = len;
    while (lo < hi) {
        const size_t mid = (lo + hi) >> 1;
        const int mid_val = ids[mid];
        if (mid_val < s)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo < len && ids[lo] == s) {
        ids.erase(ids.begin()+lo);
        vals.erase(vals.begin()+lo);
        return;
    }

}

// TODO: isValid returns true only for valid oblivious routing table, e.g. all linear flows must be stored
//   Hence, if the table should also return true for tables that only satisfy the flow conservation constraints.
bool LinearRoutingTable::isValid(const IGraph& g) const {
    const int m = src_ids.size();
    for (int e = 0; e < m; ++e) {
        const auto& ids = src_ids[e];
        for (size_t i = 1; i < ids.size(); ++i) {
            if (ids[i] <= ids[i-1]) {
                return false;
            }
        }
    }

    // check flow conservation constraint
    std::vector<double> net_flow(n, 0.0);
    for (int e = 0; e < m; ++e) {
        // get the orientation of the flow
        const auto& [u, v] = g.getEdgeEndpoints(e);
        // int sign = (u < v) ? 1 : -1;
        for (int s = 1; s < n; ++s) { // By default the root node is 0
            if ( s == 9) {
                bool stop = true;
            }
            if ( u == s) {
                double flow = getFlow(e, s);
                net_flow[s] += flow;
            }
        }
    }

    bool all_unit_flow = true;
    for (int s = 1; s < n; ++s) {
        if (std::abs(net_flow[s] - 1.0) > VERY_SOFT_EPS) {
            all_unit_flow &= false;
            std::cout << "Node " << s << " has net flow " << net_flow[s] << " (expected 1.0)\n";
        }
    }
    return all_unit_flow;
}



void LinearRoutingTable::printFlows(const IGraph& g) const {
    for (int s = 0; s < n; ++s) {
        std::cout << "Flows for source " << s << ":\n";
        for (int e = 0; e < g.getNumDirectedEdges(); ++e) {
            double flow = getFlow(e, s);
            if (std::abs(flow) > EPS) {
                auto [u, v] = g.getEdgeEndpoints(e);
                std::cout << "  Edge (" << u << ", " << v << "): " << flow << "\n";
            }
        }
    }
}

void LinearRoutingTable::printFlowsForSource(const IGraph& g, const int& s) const {
    std::cout << "Flows for source " << s << ":\n";
    for (int e = 0; e < g.getNumDirectedEdges(); ++e) {
        double flow = getFlow(e, s);
        if (std::abs(flow) > EPS) {
            auto [u, v] = g.getEdgeEndpoints(e);
            std::cout << "  Edge (" << u << ", " << v << "): " << flow << "\n";
        }
    }
}



void LinearRoutingScheme::printRoutingTable() const {
    routing_table.printFlows(g);
}

void LinearRoutingScheme::printFlowForSource(const int& s) const {
    routing_table.printFlowsForSource(g, s);
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

void LinearRoutingScheme::routeDemands(
    std::vector<double>& congestion,
    const demands& demands
) const {
    const int m = g.getNumDirectedEdges();
    const int n = g.getNumNodes();

    std::vector<double> directed_congestion(m, 0.0);

    for (int e = 0; e < m; ++e) {
        double flow = 0.0;

        for (int s = 0; s < n; ++s) {
            for (int t = 0; t < n; ++t) {
                if (s == t) {
                    continue;
                }

                const auto demandValueOpt = demands.getDemandValue(s, t);

                if (!demandValueOpt || *demandValueOpt <= 0.0) {
                    continue;
                }

                flow += *demandValueOpt * std::abs(
                    this->getFlow(e, s, t)
                );
            }
        }

        directed_congestion[e] = flow / g.getEdgeCapacity(e);
    }

    congestion.assign(m, 0.0);

    for (int e = 0; e < m; ++e) {
        const auto& [u, v] = g.getEdgeEndpoints(e);
        const int undirected_idx = (u < v ? e : g.getAntiEdge(e));

        congestion[undirected_idx] += directed_congestion[e];
    }
}

bool LinearRoutingScheme::isValid() {
    return routing_table.isValid(g);
}