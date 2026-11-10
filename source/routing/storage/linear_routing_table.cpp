//
// Created by Mert Biyikli on 23.06.26.
//
#include "routing/storage/linear_routing_table.h"
#include "utils/my_math.h"
#include <cassert>
#include <iostream>

void LinearRoutingTable::init(const optimized::Graph<EdgeData>& g) {
    const int numEdges = g.globalEdgeCount() * 2;

    n = g.getNumNodes();

    src_ids.assign(numEdges, {});
    src_flows.assign(numEdges, {});
}

void LinearRoutingTable::addFlow(int e, int s, double flow_sx) {
    assert(e >= 0 && e < src_ids.size() && s >= 0 && s < n);
    auto& ids = src_ids[e];
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
    } else {
        // insert new value at the new position
        ids.insert(ids.begin() + static_cast<long>(lo), s);
        vals.insert(vals.begin() + static_cast<long>(lo), flow_sx);
    }
}

const double LinearRoutingTable::getFlow(int e, int s) const {
    assert(e >= 0 && e < src_ids.size() && s >= 0 && s < n);
    const auto& ids = src_ids[e];
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
    auto& ids = src_ids[e];
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
        ids.erase(ids.begin() + index_s);
        vals.erase(vals.begin() + index_s);
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
        ids.erase(ids.begin() + lo);
        vals.erase(vals.begin() + lo);
        return;
    }
}

// TODO: isValid returns true only for valid oblivious routing table, e.g. all linear flows must be stored
//   Hence, if the table should also return true for tables that only satisfy the flow conservation constraints.
bool LinearRoutingTable::isValid(const optimized::Graph<EdgeData>& g) const {
    const int num_nodes = g.getNumNodes();
    const int num_edges = g.globalEdgeCount() * 2;

    if (n != num_nodes || src_ids.size() != static_cast<size_t>(num_edges) || src_flows.size() != src_ids.size()) {
        return false;
    }

    // 1. Verify the sparse routing-table structure.
    for (int e = 0; e < num_edges; ++e) {
        const auto& ids = src_ids[e];
        const auto& flows = src_flows[e];

        if (ids.size() != flows.size()) {
            return false;
        }

        for (size_t i = 0; i < ids.size(); ++i) {
            if (ids[i] < 0 || ids[i] >= n) {
                return false;
            }

            if (i > 0 && ids[i] <= ids[i - 1]) {
                return false;
            }

            if (!std::isfinite(flows[i])) {
                return false;
            }
        }
    }

    // 2. Verify the divergence of every source-to-root basis flow.
    // Root node is 0.
    for (int s = 1; s < n; ++s) {
        std::vector<double> divergence(n, 0.0);

        // Physical edge IDs are (2k, 2k+1).
        for (int e = 0; e < num_edges; e += 2) {
            if (!g.edgeAlive(e)) {
                continue;
            }

            const int reverse = g.reverse(e).id;
            const auto [u, v] = g.getEdgeEndpoints(e);

            // Signed flow in orientation u -> v.
            const double flow = getFlow(e, s) - getFlow(reverse, s);

            if (!std::isfinite(flow)) {
                return false;
            }

            divergence[u] += flow;
            divergence[v] -= flow;
        }

        std::vector<std::string> errors;

        // 3. Every vertex must satisfy the conservation constraint.
        for (int v = 0; v < n; ++v) {
            double expected = 0.0;

            if (v == s) {
                expected = 1.0;
            } else if (v == 0) {
                expected = -1.0;
            }

            if (std::abs(divergence[v] - expected) > VERY_SOFT_EPS) {
                errors.push_back("Invalid basis flow source=" + std::to_string(s) + " vertex=" + std::to_string(v) + " divergence=" + std::to_string(divergence[v]) + " expected=" + std::to_string(expected));
            }
        }
        if (errors.size()) {
            for (const auto& error : errors) {
                //std::cout << error << "\n";
            }
            return false;
        }
    }


    return true;
}

void LinearRoutingTable::printFlows(const optimized::Graph<EdgeData>& g) const {
    for (int s : g) {
        std::cout << "Flows for source " << s << ":\n";
        for (int e = 0; e < g.getNumDirectedEdges(); ++e) {
            auto [u, v] = g.getEdgeEndpoints(e);
            if (u > v) {
                continue; // print each physical edge once
            }

            const int reverse_edge = g.reverse(e).id;
            const double net_flow = getFlow(e, s) - getFlow(reverse_edge, s);

            if (std::abs(net_flow) <= EPS) {
                continue;
            }

            if (net_flow > 0.0) {
                std::cout << "  Edge (" << u << ", " << v << "): " << net_flow << "\n";
            } else {
                std::cout << "  Edge (" << v << ", " << u << "): " << -net_flow << "\n";
            }
        }
    }
}

void LinearRoutingTable::printFlowsForSource(const optimized::Graph<EdgeData>& g, const int& s) const {
    std::cout << "Flows for source " << s << ":\n";
    for (int e = 0; e < g.getNumDirectedEdges(); ++e) {
        auto [u, v] = g.getEdgeEndpoints(e);
        if (u > v) {
            continue; // print each physical edge once
        }

        const int reverse_edge = g.reverse(e).id;
        const double net_flow = getFlow(e, s) - getFlow(reverse_edge, s);

        if (std::abs(net_flow) > EPS) {
            if (net_flow > 0.0) {
                std::cout << "  Edge (" << u << ", " << v << "): " << net_flow << "\n";
            } else {
                std::cout << "  Edge (" << v << ", " << u << "): " << -net_flow << "\n";
            }
        }
    }
}

void LinearRoutingScheme::printRoutingTable() const {
    routing_table.printFlows(g);
}

void LinearRoutingScheme::printFlowForSource(const int& s) const {
    routing_table.printFlowsForSource(g, s);
}

void LinearRoutingScheme::getEdgeBasis(int edge, std::vector<double>& basis) const {

    const int reverse_edge = g.reverse(edge).id;
    const int n = g.getNumNodes();
    basis.assign(n, 0.0);

    for (int s = 0; s < n; ++s) {
        basis[s] = routing_table.getFlow(edge, s) - routing_table.getFlow(reverse_edge, s);
    }
}

double LinearRoutingScheme::computeObliviousRatio() const {
    // for the liner routing scheme, we can compute the oblivious ratio, by pushing for each
    // edge the capacity of along the edge points
    const int m = g.getNumDirectedEdges();
    demands worst_case_demands;
    for (int e = 0; e < m; ++e) {
        auto [u, v] = g.getEdgeEndpoints(e);
        if (u > v)
            continue; // only consider one orientation for undirected edges
        double capacity = g.edgeData(e).capacity;
        worst_case_demands.addDemand(u, v, capacity);
    }

    std::vector<double> congestion;
    routeDemands(congestion, worst_case_demands);
    double max_ratio = getMaxCongestion(congestion);
    return max_ratio;
}

void LinearRoutingScheme::initRoutingTable() {
    routing_table.init(g);
}

// Basis flows can be stored on either directed orientation of the same physical
// edge. Reconstruct signed basis flow on orientation e by subtracting the
// anti-edge contribution before applying linearity across sources.
double LinearRoutingScheme::getFlow(int e, int s, int t) const {
    const int reverse_edge = g.reverse(e).id;

    const double flow_sx = routing_table.getFlow(e, s) - routing_table.getFlow(reverse_edge, s);

    const double flow_tx = routing_table.getFlow(e, t) - routing_table.getFlow(reverse_edge, t);

    const double total_flow = flow_sx - flow_tx;
    return (std::abs(total_flow) < EPS ? 0.0 : total_flow);
}

double LinearRoutingScheme::getPhysicalLinkFlow(int e, int s, int t) const {
    return std::abs(getFlow(e, s, t));
}

void LinearRoutingScheme::addFlow(int e, int s, int t, double flow_sx) {
    routing_table.addFlow(e, s, flow_sx);
    routing_table.addFlow(e, t, -flow_sx);
}

void LinearRoutingScheme::routeDemands(std::vector<double>& congestion, const demands& demand_map) const {
    // Allocate according to the complete stable edge-ID space.
    const int m = g.globalEdgeCount() * 2;
    const int n = g.getNumNodes();

    congestion.assign(m, 0.0);

    std::vector<double> basis(n);

    const auto combinedDemands = combineDemandPairs(demand_map, n);

    // Physical links have directed IDs (2k, 2k+1).
    for (int e = 0; e < m; e += 2) {

        // Removed links must not carry traffic.
        if (!g.edgeAlive(e))
            continue;

        const int reverse = g.reverse(e).id;

        getEdgeBasis(e, basis);

        double load = 0.0;

        for (const auto& demand : combinedDemands) {
            const double magnitude = std::abs(basis[demand.source] - basis[demand.target]);

            if (magnitude < EPS)
                continue;

            load += demand.value * magnitude;
        }

        const double forwardCapacity = g.edgeData(e).capacity;

        const double reverseCapacity = g.edgeData(reverse).capacity;

        congestion[e] = load / forwardCapacity;
        congestion[reverse] = load / reverseCapacity;
    }
}

bool LinearRoutingScheme::isValid() {
    return routing_table.isValid(g);
}