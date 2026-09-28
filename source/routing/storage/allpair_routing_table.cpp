//
// Created by Mert Biyikli on 23.06.26.
//

#include "routing/storage/allpair_routing_table.h"
#include <cassert>
#include "utils/my_math.h"
#include "utils/hash.h"
#include <iostream>
#include <unordered_map>

void AllPairRoutingTable::init(const optimized::Graph<EdgeData>& g) {
    const int numEdges = g.getNumDirectedEdges();
    n = g.getNumNodes();
    adj_ids.assign(numEdges, {});
    adj_vals.assign(numEdges, {});
    anti_edge.resize(numEdges, INVALID_EDGE_ID);
    for (int e = 0; e < numEdges; ++e) {
        const int& anti_e = g.reverse(e).id;
        anti_edge[e] = anti_e;
    }
}


void AllPairRoutingTable::addFlow(const int& e, const int& s, const int& t, double fraction) {
    assert(e >= 0 && e < (int)adj_vals.size());
    if (s >= t) return;

    const int commodity_id = getCommodityID(n, s, t);
    if (commodity_id == INVALID_COMMODITY_ID) return;

    if (std::abs(fraction) <= SOFT_EPS) return;

    addFlow(e, commodity_id, fraction);
}



void AllPairRoutingTable::addFlow(const int e, const int commodity_id, const double delta) {
    assert(e >= 0 && e < (int)adj_vals.size());
    auto& ids  = adj_ids[e];
    auto& vals = adj_vals[e];
    const int len = (int)ids.size();

    int idx = findIndexSorted(ids, commodity_id);
    if (idx >= 0) {
        vals[idx] += delta;
        if (std::abs(vals[idx]) <= SOFT_EPS) {
            eraseAt(e, idx);
        }
        return;
    }

    // insert
    int lo = 0, hi = len;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (ids[mid] < commodity_id) lo = mid + 1;
        else hi = mid;
    }
    ids.insert(ids.begin() + lo, commodity_id);
    vals.insert(vals.begin() + lo, delta);
    if (std::abs(delta) <= SOFT_EPS) { // don't keep near-zero
        eraseAt(e, lo);
    }
}

void AllPairRoutingTable::eraseAt(int e, int idx) {
    assert(e >= 0 && e < (int)adj_vals.size());
    auto& ids  = adj_ids[e];
    auto& vals = adj_vals[e];
    ids.erase(ids.begin() + idx);
    vals.erase(vals.begin() + idx);
}

int AllPairRoutingTable::findIndexSorted(const std::vector<int>& ids, int commodity_id) const {
    const int len = (int)ids.size();
    int linear_bound = 8;
    for (int i = 0; i < std::min(len, linear_bound); ++i)
        if (ids[i] == commodity_id) return i;

    int lo = 0, hi = len;
    while (lo < hi) {
        int mid = (lo + hi) >> 1;
        if (ids[mid] < commodity_id) lo = mid + 1;
        else hi = mid;
    }
    return (lo < len && ids[lo] == commodity_id) ? lo : -1;
}




double AllPairRoutingTable::getFlow(int e , int s, int t) const {
    if ( s== t ) return 0.0;
    if (s > t) {
        return (-getFlow(e, t, s));
    }
    assert(e >= 0 && e < adj_vals.size());
    auto& ids  = adj_ids[e];
    auto& vals = adj_vals[e];
    int len = static_cast<int>(ids.size());


    int commodity_id = getCommodityID(n, s, t);
    assert(commodity_id != INVALID_COMMODITY_ID);
    // first run linear scan
    int linear_bound = 8;
    for (int i = 0; i < std::min(len, linear_bound); ++i) {
        if (ids[i] == commodity_id) {
            return vals[i];
        }
    }

    int lo = 0, hi = len;
    while (lo < hi) {
        const size_t mid = (lo + hi) >> 1;
        const auto& mid_val = ids[mid];
        if (mid_val < commodity_id)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo < len && ids[lo] == commodity_id) {
        return vals[lo];
    } else {
        return 0.0;
    }
}



bool AllPairRoutingTable::isValid(const optimized::Graph<EdgeData>& g) const {
    const int m = adj_ids.size();
    for (int e = 0; e < m; ++e) {
        const auto& ids = adj_ids[e];
        for (size_t i = 1; i < ids.size(); ++i) {
            if (ids[i] <= ids[i-1]) {
                return false;
            }
        }
    }

    // check flow conservation constraint
    const int n = g.getNumNodes();
    std::unordered_map<std::pair<int, int>, double, PairHash> net_flow;
    for (int e = 0; e < m; ++e) {
        // get the orientation of the flow
        const auto& [u, v] = g.getEdgeEndpoints(e);

        for (int s = 0; s < n; ++s) { // By default the root node is 0
            for ( int t = 0; t < n; ++t ) {
                if ( s == t ) continue;

                if ( u == s ) {
                    double flow = getFlow(e, s, t);
                    net_flow[{s, t}] += flow;
                }
                if ( v == s ) {
                    double flow = getFlow(e, s, t);
                    net_flow[{s, t}] += flow;
                }
            }
        }
    }

    bool all_unit_flow = true;
    for (int s = 0; s < n; ++s) {
        for (int t = 0; t < n; ++t) {
            if ( s == t ) continue;

            if (s < t) {
                if (std::abs(net_flow[{s, t}] - 1) > SOFT_EPS) {
                    all_unit_flow &= false;
                    std::cout << "Commodity " << s << " -> " << t << " has net flow "
                              << net_flow[{s, t}] << ")\n";
                }
            }else {
                if (std::abs(net_flow[{t, s}] - 1) > SOFT_EPS) {
                    all_unit_flow &= false;
                    std::cout << "Commodity " << t << " -> " << s << " has net flow "
                              << net_flow[{t, s}] << ")\n";
                }
            }
        }
    }
    return all_unit_flow;
}

void AllPairRoutingTable::printFlows(const optimized::Graph<EdgeData>& g) const {
    for (int s = 0; s < g.getNumNodes(); ++s) {
        for (int t = 0; t < g.getNumNodes(); ++t) {
            if (s >= t) continue;

            std::cout << "Flows for commodity " << s << " -> " << t << ":\n";
            for (int e = 0; e < g.getNumDirectedEdges(); ++e) {
                double flow = getFlow(e, s, t);
                if (std::abs(flow) > EPS) {
                    auto [u, v] = g.getEdgeEndpoints(e);
                    std::cout << "  Edge (" << u << ", " << v << "): " << flow << "\n";
                }
            }
        }
    }
}

const int AllPairRoutingTable::getSize() const {
    assert(adj_ids.size() == adj_vals.size());
    return adj_ids.size();
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

void AllPairRoutingScheme::routeDemands(
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
                    routing_table.getFlow(e, s, t)
                );
            }
        }

        directed_congestion[e] = flow / g.edgeData(e).capacity;
    }

    congestion.assign(m, 0.0);

    for (int e = 0; e < m; ++e) {
        const auto& [u, v] = g.getEdgeEndpoints(e);
        const int undirected_idx = (u < v ? e : g.reverse(e).id);

        congestion[undirected_idx] += directed_congestion[e];
    }
}


bool AllPairRoutingScheme::isValid() {
    return routing_table.isValid(g);
}
