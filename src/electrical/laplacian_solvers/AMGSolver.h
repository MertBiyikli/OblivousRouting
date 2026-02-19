//
// Created by Mert Biyikli on 06.08.25.
//

#ifndef OBLIVIOUSROUTING_AMGSOLVER_H
#define OBLIVIOUSROUTING_AMGSOLVER_H

#include <vector>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include "../../solver/routing_table.h"
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include "LaplacianSolver.h"
#include "amgcl/coarsening/runtime.hpp"
#include "amgcl/solver/runtime.hpp"
#include <amgcl/make_solver.hpp>


class AMGSolver : public LaplacianSolver {

    // we keep the components of the AMG solver as wrapper types to allow runtime configuration of coarsening and relaxation methods
    using Backend = amgcl::backend::builtin<double>;

    typedef amgcl::make_solver<
    amgcl::amg<
        Backend,
        amgcl::runtime::coarsening::wrapper,
        amgcl::runtime::relaxation::wrapper
        >,
    amgcl::runtime::solver::wrapper<Backend>
    > AMG;


    std::unique_ptr<AMG> hierarchy;
    bool use_dirichlet = true;
    int dirichlet_root = 0;
    std::vector<double> m_values_dirichlet;

public:

    std::vector<double> solve(const std::vector<double> &b, double eps = EPS) override;
    Eigen::VectorXd solve(const Eigen::VectorXd &b, double eps = EPS) override;

    void updateSolver() override;
    void buildLaplacian() override;

};

#endif //OBLIVIOUSROUTING_AMGSOLVER_H
