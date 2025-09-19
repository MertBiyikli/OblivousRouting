//
// Created by Mert Biyikli on 18.09.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_MST_SOLVER_H
#define OBLIVIOUSROUTING_RAECKE_MST_SOLVER_H

#include "raecke_random_mst.h"
#include "raecke_mst_transform.h"
#include "../raecke_framework.h"


/* * This class implements Räcke’s oblivious‐routing solver
 * as in his STOC 2008 paper using the uniform MST algorithm.
 */
using RaeckeMSTSolver = RaeckeFramework<RaeckeMST, RaeckeMSTTransform>;


#endif //OBLIVIOUSROUTING_RAECKE_MST_SOLVER_H