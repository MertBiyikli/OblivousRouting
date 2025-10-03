//
// Created by Mert Biyikli on 29.09.25.
//

#ifndef OBLIVIOUSROUTING_AMGBICGSTAB_H
#define OBLIVIOUSROUTING_AMGBICGSTAB_H

#include "LaplacianSolver.h"
#include "../../utils/hash.h"
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <memory>

#include "LaplacianSolver.h"

class AMGBiCGStabSolver : public LaplacianSolver {
private:
    using Solver = amgcl::make_solver<
        amgcl::amg<
            amgcl::backend::builtin<double>,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai0
        >,
        amgcl::solver::bicgstab<amgcl::backend::builtin<double>>
    >;


    std::unique_ptr<Solver> solver;

public:

    std::vector<double> solve(const std::vector<double>& b) override;
    Eigen::VectorXd solve(const Eigen::VectorXd& b) override;
    void updateSolver() override;
    void buildLaplacian() override;
    void buildLaplacian_(const Graph& g) override;


};

#endif //OBLIVIOUSROUTING_AMGBICGSTAB_H