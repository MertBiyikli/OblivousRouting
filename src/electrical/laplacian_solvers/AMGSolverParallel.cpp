#include "AMGSolverParallel.h"
#include <iostream>
#include <omp.h>
#include <amgcl/backend/builtin_hybrid.hpp>

#include "amgcl/solver/precond_side.hpp"


AMGSolverMT::AMGSolverMT() {
    solver_params_.solver.tol = 1e-8;
    solver_params_.solver.maxiter = 1000;

    precond_params_.npre = 1;
    precond_params_.npost = 1;
    precond_params_.ncycle = 1;
    precond_params_.pre_cycles = 1;

}

void AMGSolverMT::buildLaplacian() {
    if (debug) std::cout << "Number of threads for AMG solver pool: " << num_threads_ << "\n";

    // --- Step 1: Aggregate edge contributions into a Laplacian map ---
    std::unordered_map<std::pair<int,int>, double> L;

    for (const auto &kv : m_edge_weights) {
        int u = kv.first.first;
        int v = kv.first.second;
        double w = kv.second;
        if (u == v) continue;

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
    m_indexMap.clear();

    for (auto &entry : L) {
        int r = entry.first.first;
        int c = entry.first.second;
        double val = entry.second;
        int pos = m_row_ptr[r] + offset[r]++;
        m_col_ind[pos] = c;
        m_values[pos] = val;
        m_indexMap[{r, c}] = pos;
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
            m_indexMap[{r, row[k].first}] = start + k;
        }
    }

    // --- Step 5: Check diagonal entries ---
    if(debug) {
        for (int i = 0; i < n; ++i) {
            auto it = m_indexMap.find({i, i});
            if (it == m_indexMap.end()) {
                // If diagonal entry is missing, insert one
                std::cerr << "Warning: row " << i << " missing diagonal, inserting one.\n";
                // Insert at end of row (resizing is simpler than shifting for now)
                int insert_pos = m_row_ptr[i + 1];
                m_col_ind.insert(m_col_ind.begin() + insert_pos, i);
                m_values.insert(m_values.begin() + insert_pos, 1.0);
                for (int r = i + 1; r <= n; ++r) m_row_ptr[r]++;
                m_indexMap[{i, i}] = insert_pos;
            } else if (m_values[it->second] == 0.0) {
                std::cerr << "Warning: row " << i << " diagonal is zero, fixing to 1.\n";
                m_values[it->second] = 1.0;
            }
        }
    }


    // --- Build solver pool ---
    A_ = amgcl::backend::crs<double>(
        n, n, m_row_ptr, m_col_ind, m_values
    );

    precond_ = std::make_shared<Precond>(A_, precond_params_);
    pool_.clear();
    pool_.reserve(num_threads_);
    for (int i = 0; i < num_threads_; ++i) {
        pool_.emplace_back(precond_);
    }

    n_ = n;           // you already use 'n' in this function
    since_last_rebuild_ = 0;

    if (debug_)
        std::cout << "[AMGSolverMT] Built CRS with " << n
                  << " nodes and " << nnz << " non-zeros, initialized "
                  << num_threads_ << " solvers.\n";

}

std::vector<double> AMGSolverMT::solve(const std::vector<double>& b) {
    std::vector<double> x(b.size(), 0.0);
#ifdef _OPENMP
    int tid = omp_get_thread_num() % pool_.size();
    pool_[tid].apply(b, x);

#else
    pool_[0].apply(b, x);
#endif
    return x;
}

Eigen::VectorXd AMGSolverMT::solve(const Eigen::VectorXd& b) {
    std::vector<double> bv(b.data(), b.data() + b.size());
    auto xv = solve(bv);
    Eigen::VectorXd x(bv.size());
    for (int i = 0; i < b.size(); ++i) x[i] = xv[i];
    return x;
}

Eigen::MatrixXd AMGSolverMT::solve_many(const Eigen::MatrixXd& B) {
    const int n = B.rows();
    const int k = B.cols();
    Eigen::MatrixXd X(n, k);

    #pragma omp parallel for schedule(static)
    for (int j = 0; j < k; ++j) {
        Eigen::VectorXd rhs = B.col(j);
        X.col(j) = solve(rhs);
    }
    return X;
}



std::vector<double> AMGSolverMT::solve_single_thread(const std::vector<double>& b) {
    std::vector<double> x(b.size(), 0.0);
#ifdef _OPENMP
    int tid = omp_get_thread_num() % pool_.size();
    // ensure that only a single thread uses this solver at a time
    #pragma omp single
    {
        pool_[tid].apply(b, x);
    }
#else
    pool_[0].apply(b, x);
#endif
    return x;
}

Eigen::VectorXd AMGSolverMT::solve_single_thread(const Eigen::VectorXd& b) {
    std::vector<double> bv(b.data(), b.data() + b.size());
    auto xv = solve(bv);
    Eigen::VectorXd x(bv.size());
    for (int i = 0; i < b.size(); ++i) x[i] = xv[i];
    return x;
}

void AMGSolverMT::updateSolver() {
    A_ = amgcl::backend::crs<double>(n_, n_, m_row_ptr, m_col_ind, m_values);
    if (!precond_) {
        precond_ = std::make_shared<Precond>(A_, precond_params_);
    } else {
        auto &dst = precond_->system_matrix();
        if (A_.nrows != dst.nrows || A_.ncols != dst.ncols)
            throw std::runtime_error("Matrix structure mismatch in updateSolver()");
        std::memcpy(dst.val, A_.val, sizeof(double) * A_.nnz);
    }

    // Refresh thread clones (reuse hierarchy)
    pool_.clear();
    for (int i = 0; i < num_threads_; ++i)
        pool_.emplace_back(precond_);

    /*
    // -------------------------------------------------------------------------
    // Step 1: Bind the updated CRS view (same topology, new values)
    // -------------------------------------------------------------------------
    A_ = amgcl::backend::crs<double>(
        n_, n_, m_row_ptr, m_col_ind, m_values
    );

    // -------------------------------------------------------------------------
    // Step 2: Update or build the shared reference preconditioner once
    // -------------------------------------------------------------------------
    if (!precond_) {
        // First-time build: construct hierarchy
        precond_ = std::make_shared<Precond>(A_, precond_params_);
        if (debug_)
            std::cout << "[AMGSolverMT] Built initial preconditioner hierarchy.\n";
    } else {
        // Fastest possible value update — O(nnz)
        auto &dst = precond_->system_matrix();
        if (A_.nrows != dst.nrows || A_.ncols != dst.ncols)
            throw std::runtime_error("Matrix structure mismatch during value refresh.");

        std::memcpy(dst.val, A_.val, sizeof(double) * A_.nnz);
    }

    // -------------------------------------------------------------------------
    // Step 3: Ensure each solver in the pool has its own hierarchy
    //         and update values in-place (thread-safe)
    // -------------------------------------------------------------------------
    if (pool_.empty()) {
        pool_.resize(num_threads_);
    }

    for (int i = 0; i < static_cast<int>(pool_.size()); ++i) {
        if (!pool_[i]) {
            // Build a new solver for this thread
            pool_[i] = std::make_unique<Solver>(A_, solver_params_);
        } else {
            // Copy updated matrix values into the solver’s own backend CRS
            auto &dst = pool_[i]->system_matrix();
            if (A_.nrows != dst.nrows || A_.ncols != dst.ncols)
                throw std::runtime_error("Matrix structure mismatch during pool update.");

            std::memcpy(dst.val, A_.val, sizeof(double) * A_.nnz);
        }
    }

    // -------------------------------------------------------------------------
    // Step 4: (Optional) Debug logging
    // -------------------------------------------------------------------------
    if (debug_)
        std::cout << "[AMGSolverMT] Updated " << pool_.size()
                  << " solver(s) with new Laplacian edge weights ("
                  << A_.nnz << " nonzeros).\n";*/
}
