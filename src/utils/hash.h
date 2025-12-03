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


struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const noexcept {
        return ((size_t)p.first << 32) ^ (size_t)p.second;
    }
};

// Key types:
using TerminalPair   = std::pair<int,int>;
using DemandMap      = std::unordered_map<TerminalPair,double, PairHash>;
using EdgeDemandMap  = std::unordered_map<std::pair<int, int>, DemandMap, PairHash>;


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

struct Demand {
    int source;
    int target;

    bool operator<(const Demand& other) const {
        return std::tie(source, target) < std::tie(other.source, other.target);
    }
    bool operator==(const Demand& other) const {
        return source == other.source && target == other.target;
    }
};


namespace std {
    template<>
    struct hash<Demand> {
        std::size_t operator()(const Demand &d) const {
            return std::hash<int>()(d.source) ^ (std::hash<int>()(d.target) << 1);
        }
    };
}




#endif //OBLIVIOUSROUTING_HASH_H
