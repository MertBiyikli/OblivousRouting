//
// Created by Mert Biyikli on 12.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_MWU_RANDOM_H
#define OBLIVIOUSROUTING_RAECKE_MWU_RANDOM_H

#include "../raecke_mwu.h"
#include "raecke_oracle_random.h"

class RaeckeMWU_Random : public RaeckeMWU {
public:
    RaeckeMWU_Random(IGraph& g, int root)
        : RaeckeMWU(g, root, std::make_unique<RandomOracle>(g)) {} // TODO: replace nullptr with actual Random Oracle
};

#endif //OBLIVIOUSROUTING_RAECKE_MWU_RANDOM_H