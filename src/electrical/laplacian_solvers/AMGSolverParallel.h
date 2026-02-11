

#pragma once
#include "LaplacianSolver.h"
#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <mutex>
#include <thread>
#include <vector>
#include <amgcl/amg.hpp>
#include <Eigen/Dense>
#include <memory>
#include "../../solver/routing_table.h"
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>


template <class AMG>
struct AMGThreadClone {
    using value_type = double;
    using vector = std::vector<value_type>;

    std::shared_ptr<const AMG> shared;

    AMGThreadClone() = default;
    explicit AMGThreadClone(const std::shared_ptr<const AMG>& base)
        : shared(base)
    {}

    // Only std::vector interface for AMGCL
    void apply(const std::vector<double> &rhs, std::vector<double> &x) const {
        shared->apply(rhs, x);
    }
};


class AMGSolverMT final : public LaplacianSolver {
    using Backend = amgcl::backend::builtin<double>;

    using Precond = amgcl::amg<
        Backend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
    >;

    using Iter = amgcl::solver::cg<Backend>;

    // This is the thing that has operator()(rhs, x) and precond()
    using Solver = amgcl::make_solver<Precond, Iter>;
    std::vector<std::unique_ptr<Solver>> solvers_;

public:

    // --- Members ---
    amgcl::backend::crs<double> A_;

#ifdef _OPENMP
    int num_threads_ = omp_get_num_threads();
#else
    int num_threads_ = 1;
#endif
    bool debug_ = false;
    int n_ = 0;


    int rebuild_interval_ = 10;         // only rebuild every k calls if values change
    int since_last_rebuild_ = 0;





    AMGSolverMT();
    ~AMGSolverMT() override = default;

    Eigen::MatrixXd solve_many(const Eigen::MatrixXd& B);

    void updateSolver() override;

    // Satisfy LaplacianSolver’s abstract API
    void buildLaplacian() override;

    void set_num_threads(int t) { num_threads_ = t; }






    virtual std::vector<double> solve(const std::vector<double>& b, double eps) override {return {};};
    virtual Eigen::VectorXd solve(const Eigen::VectorXd& b, double eps) override {return {};};

};
