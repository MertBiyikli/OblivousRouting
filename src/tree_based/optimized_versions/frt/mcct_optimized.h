//
// Created by Mert Biyikli on 06.11.25.
//

#ifndef OBLIVIOUSROUTING_MCCT_OPTIMIZED_H
#define OBLIVIOUSROUTING_MCCT_OPTIMIZED_H


#include "../../../datastructures/GraphCSR.h"
#include "../../frt/mcct/Tree.h"

#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include <memory>
#include <random>
#include <cmath>
#include <iostream>

class MCCT_optimized {
public:
    struct Demand { int s, t; double w; };

    MCCT_optimized();

    void setGraph(const GraphCSR& g);
    void setSeed(uint64_t seed);
    void setMaxPairsForUnitDemands(int k);        // -1 => all s<t
    void clearDemands();                          // optional
    void addDemand(int s, int t, double w=1.0);   // optional, if you don’t want unit all-pairs

    // Build one HST using CKR/FRT procedure (1 permutation + 1 alpha).
    // Returns a valid FRT_Tree (hierarchical clustering tree).
    FRT_Tree build(bool debug=false);

private:
    GraphCSR G;
    std::mt19937_64 rng{0xC0FFEE};
    int max_pairs = -1;            // cap on #unit demands; -1 == all pairs
    std::vector<Demand> demands;

    // Helpers
    double estimateDiameter() const;                   // 2 BFS/Dijkstra heuristic
    std::vector<int> randomPermutation(int n);         // π
    double randomAlpha();                              // α ∈ [1,2)
    void truncatedDijkstra(int src, double radius,
                           std::vector<double>& dist,
                           std::vector<int>& owner,
                           int ownerId) const;

    // Construct clusters at one level given r_ℓ and permutation
    // owner[v] is assigned cluster id, cluster_id -> vector of members
    int makeLevel(double radius,
                  const std::vector<int>& perm,
                  std::vector<int>& owner,
                  std::vector<std::vector<int>>& clusters) const;

    // Build unit all-pairs if demands empty
    void ensureDemands();

    // Link clusters into HST structure
    void buildTreeFromLevels(const std::vector<std::vector<std::vector<int>>>& levels,
                             FRT_Tree& T) const;
};


#endif //OBLIVIOUSROUTING_MCCT_OPTIMIZED_H