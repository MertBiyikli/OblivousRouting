

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
public:

    // --- Backend + Templates ---
    using Backend  = amgcl::backend::builtin<double>;
    using Precond = amgcl::amg<
        Backend,
        amgcl::coarsening::smoothed_aggregation,
        amgcl::relaxation::spai0
    >;
    using Krylov   = amgcl::solver::cg<Backend>;
    using Solver   = amgcl::make_solver<Precond, Krylov>;

    // --- Members ---
    amgcl::backend::crs<double> A_;
    std::shared_ptr<Precond> precond_;
    std::vector<AMGThreadClone<Precond>> pool_;
    Precond::params precond_params_;
    Solver::params  solver_params_;
    std::mutex pool_mutex_;

    int num_threads_ = std::thread::hardware_concurrency();
    bool debug_ = false;
    int n_ = 0;


    int rebuild_interval_ = 10;         // only rebuild every k calls if values change
    int since_last_rebuild_ = 0;





    AMGSolverMT();
    ~AMGSolverMT() override = default;


    std::vector<double> solve(const std::vector<double>& b) override;  // required
    Eigen::VectorXd solve(const Eigen::VectorXd& b) override;          // convenience
    Eigen::MatrixXd solve_many(const Eigen::MatrixXd& B);

    void updateSolver() override;

    // Satisfy LaplacianSolver’s abstract API
    void buildLaplacian() override;

    void set_num_threads(int t) { num_threads_ = t; }

    void updateAllEdges(const std::vector<double>& new_weights,
                    const std::vector<std::pair<int,int>>& edges) override
    {
        if (new_weights.size() != edges.size()) {
            throw std::runtime_error("[AMGSolverMT::updateAllEdges] size mismatch between weights and edges");
        }

        // --- 1. Update Laplacian values in place ---
        for (size_t e = 0; e < edges.size(); ++e) {
            const int u = edges[e].first;
            const int v = edges[e].second;

            const double old_w = m_edge_weights[{u, v}];
            const double new_w = new_weights[e];
            const double delta = new_w - old_w;

            if (std::abs(delta) < 1e-15) continue; // skip negligible change

            // Keep symmetry
            m_edge_weights[{u, v}] = new_w;
            m_edge_weights[{v, u}] = new_w;

            // --- Laplacian updates ---
            // L(u,u) and L(v,v) increase by +delta
            if (auto it = m_indexMap.find({u, u}); it != m_indexMap.end())
                m_values[it->second] += delta;
            if (auto it = m_indexMap.find({v, v}); it != m_indexMap.end())
                m_values[it->second] += delta;

            // L(u,v) and L(v,u) decrease by -delta
            if (auto it = m_indexMap.find({u, v}); it != m_indexMap.end())
                m_values[it->second] -= delta;
            if (auto it = m_indexMap.find({v, u}); it != m_indexMap.end())
                m_values[it->second] -= delta;
        }
    }

    void buildLaplacian_(const Graph&) override {};

    std::vector<Eigen::VectorXd> solve_multi(const std::vector<Eigen::VectorXd>& rhs_batch) {
        if (rhs_batch.empty()) return {};

        const int n = rhs_batch[0].size();
        const int k = rhs_batch.size();

        std::vector<Eigen::VectorXd> results;
        results.reserve(k);

        // Reuse same solver (and hierarchy) for all RHS
        for (int j = 0; j < k; ++j) {
            Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
            std::vector<double> b(rhs_batch[j].data(), rhs_batch[j].data() + n);
            std::vector<double> x_out(n, 0.0);

            // Solve with same preconditioner
#ifdef _OPENMP
            int tid = omp_get_thread_num() % pool_.size();
            pool_[tid].apply(b, x_out);
#else
            if (pool_.empty()) {
                throw std::runtime_error("AMGSolverMT pool is empty.");
            }
            auto& solver = pool_[0];
#endif

            for (int i = 0; i < n; ++i)
                x[i] = x_out[i];

            results.push_back(std::move(x));
        }

        return results;
    }

    std::vector<double> solve_single_thread(const std::vector<double>& b);  // required
    Eigen::VectorXd solve_single_thread(const Eigen::VectorXd& b);          // convenience

};
