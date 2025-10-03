//
// Created by Mert Biyikli on 18.09.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_RANDOM_MST_H
#define OBLIVIOUSROUTING_RAECKE_RANDOM_MST_H

/**
 * This class implements Räcke’s oblivious‐routing solver
 * as in his STOC 2008 paper using the uniform MST algorithm.
 */

#include <memory>
#include <vector>
#include <unordered_map>
#include <random>
#include "mst.h"
#include "../../solver/solver.h"
#include "../raecke_base.h"

class RaeckeMST : public RaeckeBase<MSTTree> {


    // add here the random number generator
    std::random_device rd;
    std::mt19937_64 rng{rd()};
    std::uniform_int_distribution<uint64_t> dist;
    RandomMST mst;
    std::vector<int> mst_parents;

public:


    // all pure virtual function to be implemented
    void init(Graph& _g);
    double iterate(int treeIndex);
    virtual MSTTree getTree(Graph& g);
    void computeRLoads(int treeIndex,
                       MSTTree& _t,
                       Graph& copyGraph);

};

#endif //OBLIVIOUSROUTING_RAECKE_RANDOM_MST_H