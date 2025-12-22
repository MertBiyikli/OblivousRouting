//
// Created by Mert Biyikli on 08.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_MWU_CKR_H
#define OBLIVIOUSROUTING_RAECKE_MWU_CKR_H

#include "../raecke_mwu.h"
#include "raecke_oracle_ckr.h"

class RaeckeMWU_CKR : public RaeckeMWU {
public:
    RaeckeMWU_CKR(IGraph& g)
        : RaeckeMWU(g,  std::make_unique<CKROracle>(g)) {}


};


#endif //OBLIVIOUSROUTING_RAECKE_MWU_CKR_H