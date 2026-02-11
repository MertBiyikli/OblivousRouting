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
    double transformation_time = 0;

    void printTimeStats() {
        std::cout << "Solve time: " << this->solve_time << " ms\n";
        std::cout << "Transformation time: " << transformation_time << " ms\n";
        std::cout << "MWU iterations: " << this->iteration_count << "\n";
        double average_oracle_time = 0.0;
        for (double t : this->oracle_running_times) {
            average_oracle_time += t;
        }
        std::cout << "Average oracle time: " << (average_oracle_time/static_cast<double>(this->oracle_running_times.size())) << " ms\n";
    }
};

#endif //OBLIVIOUSROUTING_MWU_FRAMEWORK_H