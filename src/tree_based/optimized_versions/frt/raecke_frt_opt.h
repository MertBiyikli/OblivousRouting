//
// Created by Mert Biyikli on 05.11.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_FRT_OPT_H
#define OBLIVIOUSROUTING_RAECKE_FRT_OPT_H

#include "../raecke_optimized_base.h"
#include "raecke_frt_transform_opt.h"
// #include "mcct_optimized.h"

class RaeckeFRTOptimized : public RaeckeOptimizedBase<FRT_Tree, RaeckeFRTTransformOptimized<FRT_Tree>> {
public:
    MCCT_optimized mcct;

public:



    FRT_Tree getTree(Graph_csr& g) override;
    void computeRLoads(int treeIndex,
                       FRT_Tree& _t,
                       Graph_csr& copyGraph) override;

    void setRequirements(const Graph_csr& g);


};

#endif //OBLIVIOUSROUTING_RAECKE_FRT_OPT_H