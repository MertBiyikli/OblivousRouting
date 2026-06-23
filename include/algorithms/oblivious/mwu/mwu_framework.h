//
// Created by Mert Biyikli on 13.12.25.
//

#ifndef OBLIVIOUSROUTING_MWU_FRAMEWORK_H
#define OBLIVIOUSROUTING_MWU_FRAMEWORK_H


#include "../oblivious_solver.h"
#include "core/utils.h"


class MWUFramework : public ILinearObliviousSolverBase {
public:

    MWUFramework(IGraph& g, int root)
        : ILinearObliviousSolverBase(g, root) {}


    const int getIterationCount() const {
        return metrics.getIterationCount();
    }

    virtual void updateDistances(const std::vector<double>& distances) = 0;

    virtual void printAdditionalStats() = 0;

protected:
    MWUMetrics metrics;
};

#endif //OBLIVIOUSROUTING_MWU_FRAMEWORK_H