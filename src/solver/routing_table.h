//
// Created by Mert Biyikli on 04.12.25.
//

#ifndef OBLIVIOUSROUTING_ROUTING_TABLE_H
#define OBLIVIOUSROUTING_ROUTING_TABLE_H
#include <cassert>
#include <vector>

#include "../datastructures/IGraph.h"
#include "../utils/hash.h"

constexpr static double EPS = 1e-16;





struct EfficientRoutingTable {
    // store the flows for each commodity
    std::vector<std::vector<std::pair<int, int>>> adj_ids; // adj_ids[e] = [s1, s2, ...] list of commodities for edge e
    std::vector<std::vector<double>> adj_vals; // adj_vals[e] = [f1, f2, ...] list of flows for edge e corresponding to adj_ids

    void init(const int numEdges) {
        adj_ids.assign(numEdges, {});
        adj_vals.assign(numEdges, {});
    }

    std::vector<double> operator[](const int e) {
        assert(e >= 0 && e < adj_vals.size());
        return adj_vals[e];
    }

    void addFlow(const int e, const int s, const int t, const double fraction) {
        assert(e >= 0 && e < adj_vals.size());
        auto& ids  = adj_ids[e];
        auto& vals = adj_vals[e];
        int len = static_cast<int>(ids.size());

        // first run linear scan
        int linear_bound = 8;
        for (int i = 0; i < std::min(len, linear_bound); ++i) {
            if (ids[i] == std::make_pair(s, t)) {
                vals[i] += fraction;
                return;
            }
        }

        // binary search to find s, t in ids
        size_t lo = 0, hi = len;
        while (lo < hi) {
            const size_t mid = (lo + hi) >> 1;
            const auto& mid_val = ids[mid];
            if (mid_val < std::make_pair(s, t))
                lo = mid + 1;
            else
                hi = mid;
        }

        if (lo < len && ids[lo] == std::make_pair(s, t)) {
            vals[lo] += fraction;
        } else {
            // Usually append to the end and sort the ids w.r.t. to the terminals
            if (lo == len) {
                ids.emplace_back(s, t);
                vals.push_back(fraction);
            } else {
                ids.insert(ids.begin() + static_cast<long>(lo), {s, t});
                vals.insert(vals.begin() + static_cast<long>(lo), fraction);

            }
        }
    }

    double getFlow(int e , int s, int t) const {
        assert(e >= 0 && e < adj_vals.size());
        auto& ids  = adj_ids[e];
        auto& vals = adj_vals[e];
        int len = static_cast<int>(ids.size());


        // first run linear scan
        int linear_bound = 8;
        for (int i = 0; i < std::min(len, linear_bound); ++i) {
            if (ids[i] == std::make_pair(s, t)) {
                return vals[i];
            }
        }

        int lo = 0, hi = len;
        while (lo < hi) {
            const size_t mid = (lo + hi) >> 1;
            const auto& mid_val = ids[mid];
            if (mid_val < std::make_pair(s, t))
                lo = mid + 1;
            else
                hi = mid;
        }
        if (lo < len && ids[lo] == std::make_pair(s, t)) {
            return vals[lo];
        } else {
            return 0.0;
        }
    }

};


// For each edge e: list of (s → flow_e(s,x))
struct LinearRoutingTable {
    int n = 0; // number of nodes
    std::vector<std::vector<int>>    src_ids;   // src_ids[e]   = [s1, s2, ...]
    std::vector<std::vector<double>> src_flows; // src_flows[e] = [f_e(s1,x), f_e(s2,x), ...]

    void init(const IGraph& g) {
        const int numEdges = g.getNumEdges();
        n = g.getNumNodes();
        src_ids.assign(numEdges, {});
        src_flows.assign(numEdges, {});
    }


    void addFlow(int e, int s, double flow_sx) {
        assert(e >= 0 && e < src_ids.size() && s >= 0 && s < n);
        auto& ids  = src_ids[e];
        auto& vals = src_flows[e];

        bool found = false;
        if (e == 10) {
            found = true;
        }


        // keep the ids sorted by s for binary search later
        size_t len = ids.size();
        size_t lo = 0, hi = len;


        // 1) Linear warm-up for small lists (fast for MWU!)
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

    // return flow for unit demand s→x on edge e, or 0 if not present
    const double getFlow(int e, int s) const {
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
};



class RoutingScheme{
protected:
    const IGraph& g;
public:
    explicit RoutingScheme(const IGraph& _g):g(_g) {
    }
    virtual ~RoutingScheme() = default;

    // avoid copying
    RoutingScheme(const RoutingScheme&) = delete;
    RoutingScheme& operator=(const RoutingScheme&) = delete;

    virtual void routeDemands(std::vector<double>& congestion, const DemandMap& demands) const = 0;

    virtual double getMaxCongestion(const std::vector<double>& congestion) const {
        double max_cong = 0.0;
        for (const auto& cong : congestion) {
            if (cong > max_cong) {
                max_cong = cong;
            }
        }
        return max_cong;
    }
};


class LinearRoutingScheme : public RoutingScheme {
public:

    LinearRoutingTable routing_table;
    int root_x = 0;
    explicit LinearRoutingScheme(const IGraph& _g, int _root_x, LinearRoutingTable&& table) : RoutingScheme(_g), root_x(_root_x), routing_table(std::move(table)) {
    }

    void initRoutingTable() {
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
    double getFlow(int e, int s, int t) const {
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

        // return routing_table.getFlow(e, s) - routing_table.getFlow(e, t);
    }

    void addFlow(int e, int s, int t, double flow_sx) {
        routing_table.addFlow(e, s, flow_sx);
        routing_table.addFlow(e, t, -flow_sx);
    }

    void routeDemands(std::vector<double>& congestion,
                      const DemandMap& demands) const override
    {


        const int m = g.getNumEdges(); // undirected edges
        std::vector<double> total_flow;
        total_flow.resize(m, 0.0);

        for (int e = 0; e < m; ++e) {
            double flow = 0.0;
            // only consider sorted orientation for undirected flows
            if (g.edgeEndpoints(e).first > g.edgeEndpoints(e).second) {
                continue;
            }
            for (int d = 0; d<demands.size(); ++d) {
                const auto& [s, t] = demands.getDemandPair(d);
                const double coeff = demands.getDemandValue(d);
                if (coeff == 0.0) continue;

                flow += coeff * std::abs(getFlow(e, s, t));
            }
            total_flow[e] = flow;
        }

        // now convert the directed flows into undirected congestion
        congestion.assign(m, 0.0);
        for (int e = 0; e<m; e++) {
            if (g.edgeEndpoints(e).first < g.edgeEndpoints(e).second) {
                congestion[e] += total_flow[e];
                congestion[e] /= g.getEdgeCapacity(e);

                int anti_e = g.getAntiEdge(e);
                if (anti_e != INVALID_EDGE_ID) {
                    congestion[anti_e] = congestion[e];
                }
            }
        }
    }

};



class NonLinearRoutingScheme : public RoutingScheme {
public:
    EfficientRoutingTable routing_table;
    // to be implemented
    explicit NonLinearRoutingScheme(const IGraph& _g, EfficientRoutingTable&& table) : RoutingScheme(_g), routing_table(std::move(table)) {
    }

    double getFlow(int e, int s, int t) const{
        int e_orig = e;
        return routing_table.getFlow(e_orig, s, t);
    }

    void addFlow(int e, int s, int t, double flow_sx) {
        routing_table.addFlow(e, s,t, flow_sx);
    }

    virtual void routeDemands(std::vector<double>& congestion, const DemandMap& demands) const {
        const int m = g.getNumEdges(); // undirected edges
        std::vector<double> total_flow;
        total_flow.resize(m, 0.0);

        for (int e = 0; e < m; ++e) {
            double flow = 0.0;
            // only consider sorted orientation for undirected flows
            if (g.edgeEndpoints(e).first > g.edgeEndpoints(e).second) {
                continue;
            }
            for (int d = 0; d<demands.size(); ++d) {
                const auto& [s, t] = demands.getDemandPair(d);
                const double coeff = demands.getDemandValue(d);
                if (coeff == 0.0) continue;

                flow += coeff * std::abs(routing_table.getFlow(e, s, t));
            }
            total_flow[e] = flow;
        }

        // now convert the directed flows into undirected congestion
        congestion.assign(m, 0.0);
        for (int e = 0; e<m; e++) {
            if (g.edgeEndpoints(e).first < g.edgeEndpoints(e).second) {
                congestion[e] += total_flow[e];
                congestion[e] /= g.getEdgeCapacity(e);

                int anti_e = g.getAntiEdge(e);
                if (anti_e != INVALID_EDGE_ID) {
                    congestion[anti_e] = congestion[e];
                }
            }
        }
    }

};


#endif //OBLIVIOUSROUTING_ROUTING_TABLE_H