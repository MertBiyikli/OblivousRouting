
#include "../../include/core/routing_table.h"
#include "../../include/utils/my_math.h"

#include <cassert>
#include <iostream>
#include <map>


inline int getCommodityID(const int& n, const int& s, const int& t) {
    if (s >= t) return INVALID_COMMODITY_ID;
    return (s*n+t);
}

void AllPairRoutingTable::init(const IGraph& g) {
    const int numEdges = g.getNumDirectedEdges();
    n = g.getNumNodes();
    adj_ids.assign(numEdges, {});
    adj_vals.assign(numEdges, {});
    anti_edge.resize(numEdges, INVALID_EDGE_ID);
    for (int e = 0; e < numEdges; ++e) {
        const int& anti_e = g.getAntiEdge(e);
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

bool AllPairRoutingTable::isValid(const IGraph& g) const {
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
    std::map<std::pair<int, int>, double> net_flow;
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

void AllPairRoutingTable::printFlows(const IGraph& g) const {
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
/*
void LinearRoutingTable::reduceFlow(int e, int s, double flow) {
    assert(e >= 0 && e < src_ids.size() && s >= 0 && s < n);
    auto& ids  = src_ids[e];
    auto& vals = src_flows[e];

    // small linear scan first (your trick)
    const int len = static_cast<int>(ids.size());
    const int linear_bound = 8;
    for (int i = 0; i < std::min(len, linear_bound); ++i) {
        if (ids[i] == s) {
            vals[i] -= flow;
            if (vals[i] <= 0.0) {
                src_ids[e].erase(src_ids[e].begin() + i);
                src_flows[e].erase(src_flows[e].begin() + i);
            }
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
        vals[lo] -= flow;
        if (vals[lo] <= 0.0) {
            src_ids[e].erase(src_ids[e].begin() + lo);
            src_flows[e].erase(src_flows[e].begin() + lo);
        }
    } else {
        // flow not found, do nothing
    }
}


void LinearRoutingTable::eraseAt(int e, int s) {
    // remove flow at edge e with commodity s
    assert(e >= 0 && e < src_ids.size() && s >= 0 && s < n);
    const auto& ids  = src_ids[e];
    const auto& vals = src_flows[e];

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
        src_ids.erase(src_ids.begin()+index_s);
        src_flows.erase(src_flows.begin()+index_s);
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
        src_ids.erase(src_ids.begin()+lo);
        src_flows.erase(src_flows.begin()+lo);
        return;
    }

}

*/
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




