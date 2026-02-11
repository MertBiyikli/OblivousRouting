//
// Created by Mert Biyikli on 09.02.26.
//

#ifndef OBLIVIOUSROUTING_SOLVERTYPES_H
#define OBLIVIOUSROUTING_SOLVERTYPES_H


#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>

#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include <amgcl/solver/runtime.hpp>

#include <boost/property_tree/ptree.hpp>

using Backend = amgcl::backend::builtin<double>;

using RuntimeSolver = amgcl::make_solver<
    amgcl::amg<
        Backend,
        amgcl::runtime::coarsening::wrapper,
        amgcl::runtime::relaxation::wrapper
    >,
    amgcl::runtime::solver::wrapper<Backend>
>;


#endif //OBLIVIOUSROUTING_SOLVERTYPES_H