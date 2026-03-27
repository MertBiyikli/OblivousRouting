//
// Created by Mert Biyikli on 13.12.25.
//

#ifndef OBLIVIOUSROUTING_MWU_FRAMEWORK_H
#define OBLIVIOUSROUTING_MWU_FRAMEWORK_H

#include <iostream>
#include <vector>

class MWUFramework{
public:
    MWUFramework() {
        iteration_count = 0;
        solve_time = 0;
        transformation_time = 0;
        mwu_weight_update_time = 0;
    }

    virtual ~MWUFramework() = default;
    std::vector<double> oracle_running_times;
    int iteration_count;
    double solve_time;
    double transformation_time;
    double mwu_weight_update_time;


    virtual void updateDistances(const std::vector<double>& distances) = 0;

    void printTimeStats() {
        std::cout << "Solve time: " << this->solve_time << " micro seconds\n";
        std::cout << "Transformation time: " << transformation_time << " micro seconds\n";
        std::cout << "MWU iterations: " << this->iteration_count << "\n";
        double average_oracle_time = 0.0;
        for (double t : this->oracle_running_times) {
            average_oracle_time += t;
        }
        std::cout << "Average oracle time: " << (average_oracle_time/static_cast<double>(this->oracle_running_times.size())) << " micro seconds\n";
        std::cout << "Total MWU weight update time: " << mwu_weight_update_time << " micro seconds\n";
        printAdditionalStats();
    }

    virtual void printAdditionalStats() = 0;
};

#endif //OBLIVIOUSROUTING_MWU_FRAMEWORK_H