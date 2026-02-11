//
// Created by Mert Biyikli on 04.12.25.
//

#ifndef OBLIVIOUSROUTING_ROUTING_TABLE_H
#define OBLIVIOUSROUTING_ROUTING_TABLE_H
#include <cassert>
#include <vector>
#include <iostream>

#include "../datastructures/IGraph.h"
#include "../utils/hash.h"

constexpr static double EPS = 1e-16;
constexpr static double SOFT_EPS = 1e-3;
constexpr static double VERY_SOFT_EPS = 1e-1;

constexpr static int INVALID_COMMODITY_ID = -1;



struct RoutingTable {
    virtual ~RoutingTable() = default;
    virtual void init(const IGraph& g) = 0;

    virtual bool isValid(const IGraph& g) const = 0;
    virtual void printFlows(const IGraph& g) const = 0;
};


inline int getCommodityID(const int& n, const int& s, const int& t) {
    if (s >= t) return INVALID_COMMODITY_ID;
    return (s*n+t);
}


struct AllPairRoutingTable :public RoutingTable{

    // store the flows for each commodity
    std::vector<std::vector<int>> adj_ids; // adj_ids[e] = [s1, s2, ...] list of commodities for edge e
    std::vector<std::vector<double>> adj_vals; // adj_vals[e] = [f1, f2, ...] list of flows for edge e corresponding to adj_ids
    int n;
    std::vector<int> anti_edge; // anti_edge[e] = id of the anti-edge of e

    void init(const IGraph& g) override {
        const int numEdges = g.getNumEdges();
        n = g.getNumNodes();
        adj_ids.assign(numEdges, {});
        adj_vals.assign(numEdges, {});
        anti_edge.resize(numEdges, INVALID_EDGE_ID);
        for (int e = 0; e < numEdges; ++e) {
            const int& anti_e = g.getAntiEdge(e);
            anti_edge[e] = anti_e;
        }
    }

    const std::vector<double>& operator[](int e) const { return adj_vals[e]; }
    std::vector<double>& operator[](int e) { return adj_vals[e]; }

    void addFlow(const int& e, const int& s, const int& t, double fraction) {
        assert(e >= 0 && e < (int)adj_vals.size());
        if (s >= t) return;

        const int commodity_id = getCommodityID(n, s, t);
        if (commodity_id == INVALID_COMMODITY_ID) return;

        // cancel against anti-edge first
        cancel2Cycle(e, s, t, fraction);
        if (std::abs(fraction) <= SOFT_EPS) return;

        addFlowNo2Cycle(e, commodity_id, fraction);
    }

    bool cancel2Cycle(const int e, const int s, const int t, double& fraction) {
        const int anti_e = anti_edge[e];
        if (anti_e == INVALID_EDGE_ID) return false;

        const int commodity_id = getCommodityID(n, s, t);
        auto& idsA  = adj_ids[anti_e];
        auto& valsA = adj_vals[anti_e];

        int idx = findIndexSorted(idsA, commodity_id);
        if (idx < 0) return false;

        const double anti_flow = valsA[idx]; // since s<t stored
        if (anti_flow <= SOFT_EPS) return false;

        const double cancel = std::min(anti_flow, fraction);
        valsA[idx] -= cancel;
        fraction   -= cancel;

        if (std::abs(valsA[idx]) <= SOFT_EPS) {
            eraseAt(anti_e, idx); // erases both ids and vals
        }
        return cancel > 0.0;
    }

    void addFlowNo2Cycle(const int e, const int commodity_id, const double delta) {
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

    int findIndexSorted(const std::vector<int>& ids, int commodity_id) const {
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

    void eraseAt(int e, int idx) {
        auto& ids  = adj_ids[e];
        auto& vals = adj_vals[e];
        ids.erase(ids.begin() + idx);
        vals.erase(vals.begin() + idx);
    }





    double getFlow(int e , int s, int t) const {
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

    bool isValid(const IGraph& g) const override {
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
            const auto& [u, v] = g.edgeEndpoints(e);
            // int sign = (u < v) ? 1 : -1;
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

                if (std::abs(net_flow[{s, t}] - 1) > SOFT_EPS) {
                    all_unit_flow &= false;
                    std::cout << "Commodity " << s << " -> " << t << " has net flow "
                              << net_flow[{s, t}] << ")\n";
                }
            }

        }
        return all_unit_flow;
    }

    void printFlows(const IGraph& g) const override{

        for (int s = 0; s < g.getNumNodes(); ++s) {
            for (int t = 0; t < g.getNumNodes(); ++t) {
                if (s >= t) continue;

                std::cout << "Flows for commodity " << s << " -> " << t << ":\n";
                for (int e = 0; e < g.getNumEdges(); ++e) {
                    double flow = getFlow(e, s, t);
                    if (std::abs(flow) > EPS) {
                        auto [u, v] = g.edgeEndpoints(e);
                        std::cout << "  Edge (" << u << ", " << v << "): " << flow << "\n";
                    }
                }
            }
        }
    }
};


// For each edge e: list of (s → flow_e(s,x))
struct LinearRoutingTable : public RoutingTable {
    int n = 0; // number of nodes
    std::vector<std::vector<int>>    src_ids;   // src_ids[e]   = [s1, s2, ...]
    std::vector<std::vector<double>> src_flows; // src_flows[e] = [f_e(s1,x), f_e(s2,x), ...]

    void init(const IGraph& g) override {
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




    bool isValid(const IGraph& g) const override{
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
            const auto& [u, v] = g.edgeEndpoints(e);
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



    void printFlows(const IGraph& g) const override {
        for (int s = 0; s < n; ++s) {
            std::cout << "Flows for source " << s << ":\n";
            for (int e = 0; e < g.getNumEdges(); ++e) {
                double flow = getFlow(e, s);
                if (std::abs(flow) > EPS) {
                    auto [u, v] = g.edgeEndpoints(e);
                    std::cout << "  Edge (" << u << ", " << v << "): " << flow << "\n";
                }
            }
        }
    }
};



class RoutingScheme{
protected:
    const IGraph& g;
    RoutingTable& table;
public:
    explicit RoutingScheme(const IGraph& _g, RoutingTable& _table):g(_g), table(_table) {
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

    virtual double _getFlow(int e, int s, int t) const = 0;
    virtual void printRoutingTable() const = 0;
};


class LinearRoutingScheme : public RoutingScheme {
public:

    LinearRoutingTable routing_table;
    int root_x = 0;
    explicit LinearRoutingScheme(const IGraph& _g, int _root_x, LinearRoutingTable&& table) : RoutingScheme(_g, table), root_x(_root_x), routing_table(std::move(table)) {
    }

    void initRoutingTable() {
        routing_table.init(g);
    }

    void printRoutingTable() const override{
        routing_table.printFlows(g);
    }

    virtual double _getFlow(int e, int s, int t) const override {
        return this->getFlow(e, s, t);
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



class AllPairRoutingScheme : public RoutingScheme {
public:
    AllPairRoutingTable routing_table;
    // to be implemented
    explicit AllPairRoutingScheme(const IGraph& _g, AllPairRoutingTable&& table) : RoutingScheme(_g, table), routing_table(std::move(table)) {
    }

    virtual double _getFlow(int e, int s, int t) const override {
        return this->getFlow(e, s, t);
    }

    void printRoutingTable() const override{
        routing_table.printFlows(g);
    }

    double getFlow(int e, int s, int t) const{
        int e_orig = e;
        return routing_table.getFlow(e_orig, s, t);
    }

    void addFlow(int e, int s, int t, double flow_sx) {
        routing_table.addFlow(e, s,t, flow_sx);
    }

    virtual void routeDemands(std::vector<double>& congestion, const DemandMap& demands) const override {
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