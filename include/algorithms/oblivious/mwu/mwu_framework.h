//
// Created by Mert Biyikli on 13.12.25.
//

#ifndef OBLIVIOUSROUTING_MWU_FRAMEWORK_H
#define OBLIVIOUSROUTING_MWU_FRAMEWORK_H


#include "../oblivious_solver.h"
#include "routing/routing_result.h"


class MWUFramework : public ILinearObliviousSolverBase {
public:

    MWUFramework(optimized::Graph<EdgeData>& g, int root)
        : ILinearObliviousSolverBase(g, root) {}


    const int getIterationCount() const {
        return metrics.getIterationCount();
    }

    virtual Result<void> updateDistances(const std::vector<double>& distances) = 0;

    MWUMetrics getMetrics() const {
        return metrics;
    }
    virtual void printAdditionalStats() = 0;

protected:
    MWUMetrics metrics;
};

#endif //OBLIVIOUSROUTING_MWU_FRAMEWORK_H