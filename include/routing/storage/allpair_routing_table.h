//
// Created by Mert Biyikli on 23.06.26.
//

#ifndef OBLIVIOUSROUTING_ALLPAIR_ROUTING_TABLE_H
#define OBLIVIOUSROUTING_ALLPAIR_ROUTING_TABLE_H

#include "../routing_table.h"

struct AllPairRoutingTable :public RoutingTable{

    // store the flows for each commodity
    std::vector<std::vector<int>> adj_ids; // adj_ids[e] = [s1, s2, ...] list of commodities for edge e
    std::vector<std::vector<double>> adj_vals; // adj_vals[e] = [f1, f2, ...] list of flows for edge e corresponding to adj_ids
    std::vector<int> anti_edge; // anti_edge[e] = id of the anti-edge of e

    void init(const optimized::Graph<EdgeData>& g) override;
    const std::vector<double>& operator[](int e) const { return adj_vals[e]; }
    std::vector<double>& operator[](int e) { return adj_vals[e]; }

    void addFlow(const int& e, const int& s, const int& t, double fraction) ;
    void addFlow(const int e, const int commodity_id, const double delta);

    int findIndexSorted(const std::vector<int>& ids, int commodity_id) const;

    void eraseAt(int e, int idx);
    double getFlow(int e , int s, int t) const;

    bool isValid(const optimized::Graph<EdgeData>& g) const override;
    void printFlows(const optimized::Graph<EdgeData>& g) const override;
    const int getSize() const override;
};


class AllPairRoutingScheme : public RoutingScheme {
public:
    AllPairRoutingTable routing_table;
    // to be implemented
    explicit AllPairRoutingScheme(const optimized::Graph<EdgeData>& _g, AllPairRoutingTable&& table) : RoutingScheme(_g), routing_table(std::move(table)) {
    }

    double getFlow(int e, int s, int t) const override;

    void addFlow(int e, int s, int t, double flow_sx);

    virtual void routeDemands(std::vector<double>& congestion, const demands& demands) const override;

    virtual void printRoutingTable() const override;

    bool isValid() override;

    void printFlowForSource(const int& s) const override {};
};
#endif //OBLIVIOUSROUTING_ALLPAIR_ROUTING_TABLE_H