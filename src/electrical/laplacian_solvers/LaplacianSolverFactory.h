//
// Created by Mert Biyikli on 29.09.25.
//

#ifndef OBLIVIOUSROUTING_LAPLACIANSOLVERFACTORY_H
#define OBLIVIOUSROUTING_LAPLACIANSOLVERFACTORY_H

#include "LaplacianSolver.h"
#include "AMGSolver.h"           // AMG + CG
#include "AMGSolverParallel.h"   // AMG + Parallel CG



#include <memory>
#include <string>
#include <stdexcept>

class SolverFactory {
public:
    enum class SolverType {
        AMG_CG,
        AMG_BICGSTAB,
        AMG_PARALLEL
    };

    static std::unique_ptr<LaplacianSolver> create(SolverType type) {
        switch (type) {
            case SolverType::AMG_CG:
                return std::make_unique<AMGSolver>();
            case SolverType::AMG_PARALLEL:
                return std::make_unique<AMGSolverMT>();
            default:
                throw std::invalid_argument("Unknown solver type");
        }
    }

    static std::unique_ptr<LaplacianSolver> create(const std::string& name) {
        if (name == "amg_cg") return create(SolverType::AMG_CG);
        if (name == "amg_parallel") return create(SolverType::AMG_PARALLEL);
        if (name == "amg_bicgstab") return create(SolverType::AMG_BICGSTAB);
        throw std::invalid_argument("Unknown solver name: " + name);
    }
};


#endif //OBLIVIOUSROUTING_LAPLACIANSOLVERFACTORY_H