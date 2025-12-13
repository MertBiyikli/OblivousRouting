//
// Created by Mert Biyikli on 11.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_MWU_FRT_H
#define OBLIVIOUSROUTING_RAECKE_MWU_FRT_H
#include "../raecke_mwu.h"
#include "raecke_oracle_frt.h"

class RaeckeMWU_FRT : public RaeckeMWU {
public:
    RaeckeMWU_FRT(IGraph& g, int root)
        : RaeckeMWU(g, root, std::make_unique<FRTOracle>(g)) {}


};
#endif //OBLIVIOUSROUTING_RAECKE_MWU_FRT_H