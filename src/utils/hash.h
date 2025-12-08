//
// Created by Mert Biyikli on 10.06.25.
//

#ifndef OBLIVIOUSROUTING_HASH_H
#define OBLIVIOUSROUTING_HASH_H

#include <map>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <cassert>



struct DemandMap {
    std::vector<int> source, target;
    std::vector<double> demand_values;

    void addDemand(int s, int t, double demand) {
        assert(demand_values.size() == source.size());
        assert(demand_values.size() == target.size());
        source.push_back(s);
        target.push_back(t);
        demand_values.push_back(demand);
    }

    // Simple iterator methods over demands
    size_t size() const {
        assert(demand_values.size() == source.size());
        assert(demand_values.size() == target.size());
        return demand_values.size();
    }
    std::pair<int, int> getDemandPair(size_t idx) const {
        return {source[idx], target[idx]};
    }
    double getDemandValue(size_t idx) const {
        return demand_values[idx];
    }
};

struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const noexcept {
        return ((size_t)p.first << 32) ^ (size_t)p.second;
    }

};


// TODO: This is used in the RaeckeTransform and try to avoid using these
// Key types:
using TerminalPair   = std::pair<int,int>;
using Demand = std::pair<int, int>;
using DemandMap_      = std::unordered_map<Demand,double, PairHash>;
using EdgeDemandMap  = std::unordered_map<std::pair<int, int>, DemandMap_, PairHash>;


// adjust this to support directed edge struct: struct Arc {int src, trg, cap, operator=(
namespace std {
    template <>
    struct hash<std::tuple<int, int, int>> {
        std::size_t operator()(const std::tuple<int, int, int>& t) const {
            std::size_t h1 = std::hash<int>{}(std::get<0>(t));
            std::size_t h2 = std::hash<int>{}(std::get<1>(t));
            std::size_t h3 = std::hash<int>{}(std::get<2>(t));
            // Combine hashes
            return (((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1));
        }
    };
}

namespace std {
    template <>
    struct hash<std::tuple<int, std::pair<int, int>>> {
        std::size_t operator()(const std::tuple<int, std::pair<int, int>>& t) const {
            std::size_t h1 = std::hash<int>{}(std::get<0>(t));
            std::size_t h2 = std::hash<int>{}(std::get<1>(t).first);
            std::size_t h3 = std::hash<int>{}(std::get<1>(t).second);
            // Combine hashes
            return (((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1));
        }
    };
}








#endif //OBLIVIOUSROUTING_HASH_H
