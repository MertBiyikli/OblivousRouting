//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_SCHEME_H
#define OBLIVIOUSROUTING_ROUTING_SCHEME_H

#include "../data_structures/graph/Igraph.h"
#include "../utils/demands.h"
#include "routing_table.h"

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

    virtual void routeDemands(std::vector<double>& congestion, const demands& demands) const = 0;
    virtual double getFlow(int e, int s, int t) const = 0;
    virtual void printRoutingTable() const = 0;

    double getMaxCongestion(const std::vector<double>& congestion) const;
    virtual bool isValid() = 0;
};



class LinearRoutingScheme : public RoutingScheme {
public:

    LinearRoutingTable routing_table;
    int root_x = 0;
    explicit LinearRoutingScheme(const IGraph& _g, int _root_x, LinearRoutingTable&& table)
    : RoutingScheme(_g, table),
    root_x(_root_x),
    routing_table(std::move(table)) {
    }

    double computeObliviousRatio();

    void initRoutingTable();

    double getFlow(int e, int s, int t) const override;

    void addFlow(int e, int s, int t, double flow_sx);

    void routeDemands(std::vector<double>& congestion,
                      const demands& demands) const override;

    virtual void printRoutingTable() const override;
    bool isValid() override;
};



class AllPairRoutingScheme : public RoutingScheme {
public:
    AllPairRoutingTable routing_table;
    // to be implemented
    explicit AllPairRoutingScheme(const IGraph& _g, AllPairRoutingTable&& table) : RoutingScheme(_g, table), routing_table(std::move(table)) {
    }

    double getFlow(int e, int s, int t) const override;

    void addFlow(int e, int s, int t, double flow_sx);

    virtual void routeDemands(std::vector<double>& congestion, const demands& demands) const override;

    virtual void printRoutingTable() const override;

    bool isValid() override;
};


#endif //OBLIVIOUSROUTING_ROUTING_SCHEME_H