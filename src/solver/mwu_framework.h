//
// Created by Mert Biyikli on 13.12.25.
//

#ifndef OBLIVIOUSROUTING_MWU_FRAMEWORK_H
#define OBLIVIOUSROUTING_MWU_FRAMEWORK_H

#include "solver.h"

class MWUFramework{
public:
    std::vector<double> oracle_running_times;
    std::vector<double> pure_oracles_running_times;
    int iteration_count = 0;
    double solve_time = 0; // total time spent in solve()
};

#endif //OBLIVIOUSROUTING_MWU_FRAMEWORK_H