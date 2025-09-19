//
// Created by Mert Biyikli on 18.08.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_SOLVER_H
#define OBLIVIOUSROUTING_RAECKE_SOLVER_H

#include "../../solver/solver.h"
#include "raecke_frt.h"
#include "raecke_frt_transform.h"
#include "../raecke_base.h"
#include "../raecke_framework.h"


using RaeckeFRTSolver = RaeckeFramework<RaeckeFRT, RaeckeFRTTransform>;

#endif //OBLIVIOUSROUTING_RAECKE_SOLVER_H