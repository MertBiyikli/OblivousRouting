#include "AMGSolverParallel.h"
#include <iostream>
#include <omp.h>
#include <amgcl/backend/builtin_hybrid.hpp>
#include "../../utils/hash.h"
#include "amgcl/solver/precond_side.hpp"

/*
 * AMGSolverMT implementation
 * Multi-threaded AMG solver using a pool of threads for solving Laplacian system
 */
AMGSolverMT::AMGSolverMT() {


    #ifdef _OPENMP
        set_num_threads(omp_get_max_threads());
    #else
        set_num_threads(1);
    #endif
}


/*
 * Build the Laplacian matrix in CRS format
 *
 */
void AMGSolverMT::buildLaplacian() {
     std::cout << "Number of threads for AMG solver pool: " << num_threads_ << "\n";

    // --- Step 1: Aggregate edge contributions into a Laplacian map ---
    std::unordered_map<std::pair<int,int>, double, PairHash> L;

    for (int e = 0; e<weight_model.weights.size(); e++) {
        int u = weight_model.from[e];
        int v = weight_model.to[e];
        double w = weight_model.getEdgeWeight(u, v);
        if (u >= v) continue;

        // Laplacian contributions:
        L[{u,u}] += w;
        L[{v,v}] += w;
        L[{u,v}] -= w;
        L[{v,u}] -= w;
    }

    // --- Step 2: Count non-zeros per row ---
    m_row_ptr.assign(n + 1, 0);
    for (auto &entry : L) {
        int r = entry.first.first;
        m_row_ptr[r + 1]++;
    }
    for (int i = 0; i < n; i++) {
        m_row_ptr[i + 1] += m_row_ptr[i];
    }

    // --- Step 3: Fill col_ind and values ---
    int nnz = (int)L.size();
    m_col_ind.resize(nnz);
    m_values.resize(nnz);

    std::vector<int> offset(n, 0);


    for (auto &entry : L) {
        int r = entry.first.first;
        int c = entry.first.second;
        double val = entry.second;
        int pos = m_row_ptr[r] + offset[r]++;
        m_col_ind[pos] = c;
        m_values[pos] = val;
        weight_model.setLaplacianIndex(r, c, pos);
    }

    // --- Step 4: Sort columns within each row ---
    for (int r = 0; r < n; ++r) {
        int start = m_row_ptr[r];
        int end = m_row_ptr[r + 1];
        std::vector<std::pair<int, double>> row;
        row.reserve(end - start);
        for (int k = start; k < end; ++k)
            row.emplace_back(m_col_ind[k], m_values[k]);

        std::sort(row.begin(), row.end(), [](auto &a, auto &b) { return a.first < b.first; });

        for (int k = 0; k < (int)row.size(); ++k) {
            m_col_ind[start + k] = row[k].first;
            m_values[start + k] = row[k].second;

            weight_model.setLaplacianIndex(r, row[k].first, start + k);
        }
    }


    if (use_dirichlet) {
        m_values_dirichlet = m_values;
        applyDirichletInPlace(m_values_dirichlet);
    } else {
        m_values_dirichlet.clear();
    }
    // --- Build CRS matrix from active values ---
    const std::vector<double>& active_vals = (use_dirichlet ? m_values_dirichlet : m_values);

    A_ = amgcl::backend::crs<double>(n, n, m_row_ptr, m_col_ind, active_vals);

    solvers_.clear();
    solvers_.reserve(num_threads_);
    for (int t = 0; t < num_threads_; ++t) {
        solvers_.emplace_back(std::make_unique<Solver>(A_));
    }


    n_ = n;

}


Eigen::MatrixXd AMGSolverMT::solve_many(const Eigen::MatrixXd& B) {
    const int n = (int)B.rows();
    const int k = (int)B.cols();
    Eigen::MatrixXd X(n, k);

#ifdef _OPENMP
    omp_set_dynamic(0);
    #pragma omp single
    {
        const int tid = omp_get_thread_num();
        std::cout << "I am thread: " << tid << " processing nodes "  << "\n";

        // Thread-local buffers reused across columns
        std::vector<double> rhs(n);
        std::vector<double> sol(n);
        for (int j = 0; j < k; ++j) {
            // copy Eigen column -> rhs
            // (B.col(j) is contiguous in memory for column-major Eigen matrices)
            std::memcpy(rhs.data(), B.col(j).data(), sizeof(double) * n);

            std::fill(sol.begin(), sol.end(), 0.0);
            (*solvers_[tid])(rhs, sol);


            // copy sol -> Eigen column
            std::memcpy(X.col(j).data(), sol.data(), sizeof(double) * n);
        }
    }
#else
    std::vector<double> rhs(n), sol(n);
    for (int j = 0; j < k; ++j) {
        std::memcpy(rhs.data(), B.col(j).data(), sizeof(double) * n);
        std::fill(sol.begin(), sol.end(), 0.0);
        (*solvers_[0])(rhs, sol);
        std::memcpy(X.col(j).data(), sol.data(), sizeof(double) * n);
    }
#endif

    return X;
}




void AMGSolverMT::updateSolver() {
    // 1) Apply Dirichlet if enabled
    if (use_dirichlet) {
        m_values_dirichlet = m_values;
        applyDirichletInPlace(m_values_dirichlet);
    } else {
        m_values_dirichlet.clear();
    }

    const std::vector<double>& active_vals =
        (use_dirichlet ? m_values_dirichlet : m_values);

    // 2) Rebuild A_ view with the *same structure*, new values
    A_ = amgcl::backend::crs<double>(n_, n_, m_row_ptr, m_col_ind, active_vals);


    for (int t = 0; t < num_threads_; ++t) {
        solvers_.emplace_back(std::make_unique<Solver>(A_));
    }
}


