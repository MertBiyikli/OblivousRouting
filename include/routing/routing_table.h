//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_ROUTING_TABLE_H
#define OBLIVIOUSROUTING_ROUTING_TABLE_H

#include "../data_structures/graph/Igraph.h"
#include "../utils/demands.h"

/*
 * Helpers:
 */
constexpr static int INVALID_COMMODITY_ID = -1;

inline int getCommodityID(const int& n, const int& s, const int& t) {
    if (s >= t) return INVALID_COMMODITY_ID;
    return (s*n+t);
}


// Base classes for storing the routing tables
class RoutingTable {
public:
    int n=0;

    virtual ~RoutingTable() = default;
    virtual void init(const IGraph& g) = 0;

    virtual bool isValid(const IGraph& g) const = 0;
    virtual void printFlows(const IGraph& g) const = 0;

    const int getNumNodes() const { return n; }
    virtual const int getSize() const = 0;
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

    virtual void routeDemands(std::vector<double>& congestion, const demands& demands) const = 0;
    virtual double getFlow(int e, int s, int t) const = 0;
    virtual void printRoutingTable() const = 0;
    virtual void printFlowForSource(const int& s) const = 0;

    virtual bool isValid() = 0;

    double getMaxCongestion(const std::vector<double>& congestion) const {
        double max_cong = 0.0;
        for (const auto& cong : congestion) {
            if (cong > max_cong) {
                max_cong = cong;
            }
        }
        return max_cong;
    }
};
#endif //OBLIVIOUSROUTING_ROUTING_TABLE_H