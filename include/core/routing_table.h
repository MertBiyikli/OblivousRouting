//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_TABLE_H
#define OBLIVIOUSROUTING_ROUTING_TABLE_H

#include "../data_structures/graph/Igraph.h"

constexpr static int INVALID_COMMODITY_ID = -1;

class RoutingTable {
public:
    virtual ~RoutingTable() = default;
    virtual void init(const IGraph& g) = 0;

    virtual bool isValid(const IGraph& g) const = 0;
    virtual void printFlows(const IGraph& g) const = 0;
};

struct AllPairRoutingTable :public RoutingTable{

    // store the flows for each commodity
    std::vector<std::vector<int>> adj_ids; // adj_ids[e] = [s1, s2, ...] list of commodities for edge e
    std::vector<std::vector<double>> adj_vals; // adj_vals[e] = [f1, f2, ...] list of flows for edge e corresponding to adj_ids
    int n;
    std::vector<int> anti_edge; // anti_edge[e] = id of the anti-edge of e

    void init(const IGraph& g) override;
    const std::vector<double>& operator[](int e) const { return adj_vals[e]; }
    std::vector<double>& operator[](int e) { return adj_vals[e]; }

    void addFlow(const int& e, const int& s, const int& t, double fraction) ;
    void addFlow(const int e, const int commodity_id, const double delta);

    int findIndexSorted(const std::vector<int>& ids, int commodity_id) const;

    void eraseAt(int e, int idx);
    double getFlow(int e , int s, int t) const;

    bool isValid(const IGraph& g) const override;
    void printFlows(const IGraph& g) const override;
};

// For each edge e: list of (s → flow_e(s,x))
class LinearRoutingTable : public RoutingTable {
public:
    int n = 0; // number of nodes
    std::vector<std::vector<int>>    src_ids;   // src_ids[e]   = [s1, s2, ...]
    std::vector<std::vector<double>> src_flows; // src_flows[e] = [f_e(s1,x), f_e(s2,x), ...]

    void init(const IGraph& g) override;
    void addFlow(int e, int s, double flow_sx);
    void eraseFlow(int e, int s);

    // return flow for unit demand s→x on edge e, or 0 if not present
    const double getFlow(int e, int s) const;
    bool isValid(const IGraph& g) const override;
    void printFlows(const IGraph& g) const override;

};
#endif //OBLIVIOUSROUTING_ROUTING_TABLE_H