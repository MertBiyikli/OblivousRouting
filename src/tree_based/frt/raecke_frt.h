//
// Created by Mert Biyikli on 05.06.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_TREE_DECOMP_H
#define OBLIVIOUSROUTING_RAECKE_TREE_DECOMP_H


#include <memory>
#include <vector>
#include <unordered_map>
#include "mcct/mcct_derandomized_weighted_solver.h"
#include "../../solver/solver.h"
#include "../raecke_base.h"



/**
 * This class implements Räcke’s oblivious‐routing solver
 * as in his STOC 2008 paper using the FRT algorithm.
 */

class RaeckeFRT : public RaeckeBase<FRT_Tree> {

    // TODO: for now keep an extra vector that is managed by the RaeckeFRT class
    //       later we can move this to the Solver class



    MCCT_Solver m_mcct;


    void init(Graph& g);
public:



    FRT_Tree getTree(Graph& g);
    void computeRLoads(int treeIndex,
                       FRT_Tree& _t,
                       Graph& copyGraph);


    void setRequirements(const Graph& g);

};





#endif //OBLIVIOUSROUTING_RAECKE_TREE_DECOMP_H
