//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_HASH_H
#define OBLIVIOUSROUTING_HASH_H

struct PairHash {
    size_t operator()(const std::pair<int,int>& p) const noexcept {
        return ((size_t)p.first << 32) ^ (size_t)p.second;
    }
};

struct TuplePair {
    std::size_t operator()(const std::tuple<int, std::pair<int, int>>& t) const {
        std::size_t h1 = std::hash<int>{}(std::get<0>(t));
        std::size_t h2 = std::hash<int>{}(std::get<1>(t).first);
        std::size_t h3 = std::hash<int>{}(std::get<1>(t).second);
        // Combine hashes
        return (((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1));
    }
};

struct TrippleTuple {
    std::size_t operator()(const std::tuple<int, int, int>& t) const {
        std::size_t h1 = std::hash<int>{}(std::get<0>(t));
        std::size_t h2 = std::hash<int>{}(std::get<1>(t));
        std::size_t h3 = std::hash<int>{}(std::get<2>(t));
        // Combine hashes
        return (((h1 ^ (h2 << 1)) >> 1) ^ (h3 << 1));
    }
};

#endif //OBLIVIOUSROUTING_HASH_H