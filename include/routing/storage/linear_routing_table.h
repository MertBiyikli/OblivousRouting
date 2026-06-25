//
// Created by Mert Biyikli on 23.06.26.
//

#ifndef OBLIVIOUSROUTING_LINEAR_ROUTING_TABLE_H
#define OBLIVIOUSROUTING_LINEAR_ROUTING_TABLE_H

#include "../routing_table.h"


// For each edge e: list of (s → flow_e(s,x))
class LinearRoutingTable : public RoutingTable {
public:
    std::vector<std::vector<int>>    src_ids;   // src_ids[e]   = [s1, s2, ...]
    std::vector<std::vector<double>> src_flows; // src_flows[e] = [f_e(s1,x), f_e(s2,x), ...]

    void init(const IGraph& g) override;
    void addFlow(int e, int s, double flow_sx);
    void eraseFlow(int e, int s);

    // return flow for unit demand s→x on edge e, or 0 if not present
    const double getFlow(int e, int s) const;
    bool isValid(const IGraph& g) const override;
    void printFlows(const IGraph& g) const override;
    void printFlowsForSource(const IGraph& g, const int& source) const;
    const int getSize() const override;
};


class LinearRoutingScheme : public RoutingScheme {
public:

    LinearRoutingTable routing_table;
    int root_x = 0;
    explicit LinearRoutingScheme(const IGraph& _g, int _root_x, LinearRoutingTable&& table)
    : RoutingScheme(_g),
    root_x(_root_x),
    routing_table(std::move(table)) {
    }

    double computeObliviousRatio() const;

    void initRoutingTable();

    double getFlow(int e, int s, int t) const override;

    void addFlow(int e, int s, int t, double flow_sx);

    void routeDemands(std::vector<double>& congestion,
                      const demands& demands) const override;

    virtual void printRoutingTable() const override;

    bool isValid() override;
    virtual void printFlowForSource(const int& s) const override;
};


#endif //OBLIVIOUSROUTING_LINEAR_ROUTING_TABLE_H