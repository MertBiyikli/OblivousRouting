//
// Created by Mert Biyikli on 06.08.25.
//

#include "AMGSolver.h"
#include <algorithm>
#include "../../utils/hash.h"
#include "amgcl/solver/cg.hpp"


void AMGSolver::buildLaplacian() {
    std::unordered_map<std::pair<int,int>, double, PairHash> L;

    for (int e = 0; e<weight_model.weights.size(); e++) {
        int u = weight_model.from[e];
        int v = weight_model.to[e];
        double w = weight_model.getEdgeWeight(u, v);
        if (u >= v) continue;

        L[{u,u}] += w;
        L[{v,v}] += w;
        L[{u,v}] -= w;
        L[{v,u}] -= w;
    }

    // Ensure every row has an explicit diagonal entry (needed for isolated vertices)
    for (int i = 0; i < n; ++i) {
        // If the diagonal is missing, this creates it with value 0.0.
        // (operator[] default-constructs to 0 and then we add 0)
        L[{i,i}] += 0.0;
    }

    const double eps_diag = 1e-12;
    for (int i = 0; i < n; ++i) {
        std::unordered_map<std::pair<int, int>, double, PairHash>::iterator it = L.find({i, i});
        if (it == L.end()) {
            L[{i,i}] = eps_diag;
        } else if (it->second == 0.0) {
            // Only bump true isolates (degree 0 => diagonal 0)
            it->second = eps_diag;
        }
    }




    // Count non-zeros per row
    m_row_ptr.assign(n + 1, 0);
    for (auto &entry : L) {
        int r = entry.first.first;
        m_row_ptr[r + 1]++;
    }
    for (int i = 0; i < n; i++) {
        m_row_ptr[i + 1] += m_row_ptr[i];
    }

    // Fill col_ind and values
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


    if(debug){
        if (!checkMatrix()) {
            throw std::runtime_error("Matrix structure is invalid");
        }
    }


    assert(checkMatrix() && "Matrix structure is invalid during buildLaplacian");

    if (use_dirichlet) {
        m_values_dirichlet = m_values;
        applyDirichletInPlace(m_values_dirichlet);
        hierarchy = std::make_unique<AMG>(std::tie(n, m_row_ptr, m_col_ind, m_values_dirichlet), params);
    } else {
        hierarchy = std::make_unique<AMG>(std::tie(n, m_row_ptr, m_col_ind, m_values), params);
    }


}



std::vector<double> AMGSolver::solve(const std::vector<double> &b, double eps) {
    if (b.size() != (size_t)n)
        throw std::runtime_error("RHS size mismatch");

    if (!hierarchy) throw std::runtime_error("AMGSolver: hierarchy not built");

    x_buffer.assign(n, 0.0);
    bvec_buffer.assign(b.begin(), b.end());

    if (use_dirichlet) {
        bvec_buffer[dirichlet_root] = 0.0;
        x_buffer[dirichlet_root] = 0.0;
    } else {
        double mean_b = std::reduce(bvec_buffer.begin(), bvec_buffer.end(), 0.0) / n;
        for (double &v : bvec_buffer) v -= mean_b;
    }

    amgcl::solver::cg<Backend>::params solver_params;
    solver_params.tol = eps;
    solver_params.maxiter = 1000; // optional

    const amgcl::solver::cg<Backend> solver(n, solver_params);

    if (debug) {
        const auto res = solver(*hierarchy, bvec_buffer, x_buffer);
        std::cout << "[AMGSolver] Iterations: " << std::get<0>(res)
                  << "  Error: " << std::get<1>(res) << std::endl;
    } else {
        solver(*hierarchy, bvec_buffer, x_buffer);
    }

    if (use_dirichlet) {
        x_buffer[dirichlet_root] = 0.0;
    } else {
        double mean_x = std::reduce(x_buffer.begin(), x_buffer.end(), 0.0) / n;
        for (double &v : x_buffer) v -= mean_x;
    }

    return x_buffer;
}


Eigen::VectorXd AMGSolver::solve(const Eigen::VectorXd& b, double eps) {
    const int n = b.size();
    std::vector<double> bvec(n);
    std::memcpy(bvec.data(), b.data(), n * sizeof(double));

    result.clear();
    result = this->solve(bvec);
    Eigen::VectorXd eigen_output(result.size());
    std::memcpy(eigen_output.data(), result.data(), n * sizeof(double));
    return eigen_output;
}


void AMGSolver::updateSolver() {
    assert(checkMatrix() && "Matrix structure is invalid during updateSolver");
    // Update hierarchy numeric values (same structure)
    if (use_dirichlet) {
        m_values_dirichlet = m_values;
        applyDirichletInPlace(m_values_dirichlet);
        hierarchy->precond().rebuild(std::tie(n, m_row_ptr, m_col_ind, m_values_dirichlet));
    } else {
        hierarchy->precond().rebuild(std::tie(n, m_row_ptr, m_col_ind, m_values));
    }
}




