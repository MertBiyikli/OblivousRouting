//
// Created by Mert Biyikli on 20.03.26.
//

#include "../../../../../include/algorithms/mwu/oracle/electrical/laplacian_solver.h"
#include "../../../../../include/utils/my_math.h"
#include "../../../../../include/utils/hash.h"
#include <amgcl/backend/builtin.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include "amgcl/solver/runtime.hpp"
#include <amgcl/make_solver.hpp>
#include <cassert>


void LaplacianSolver::init(IGraph& g, std::vector<double>& _adj_edge_weights, int n, const std::vector<std::pair<int, int>>& edges, bool debug) {
     assert(g.getNumUndirectedEdges() == static_cast<int>(_adj_edge_weights.size()));
     this->debug = debug;
     this->n = n;
     m_row_ptr.clear();
     m_col_ind.clear();
     m_values.clear();

     weight_model.init(g);
     for (int e = 0; e < edges.size(); e++) {
         int u = edges[e].first;
         int v = edges[e].second;
         double w = _adj_edge_weights[e];
         weight_model.setEdgeWeight(u, v, w);
         weight_model.setEdgeWeight(v, u, w); // keep it symmetric
     }

     buildLaplacian();
 }

void LaplacianSolver::updateAllEdges(const std::vector<double> &new_weights, const std::vector<std::pair<int, int> > &edges) {
     assert(new_weights.size() == edges.size());

     for (size_t e = 0; e < edges.size(); ++e) {
         int u = edges[e].first;
         int v = edges[e].second;
         double old_w = weight_model.getEdgeWeight(u, v);
         double new_w = std::max(new_weights[e], EPS);
         double delta = new_w - old_w;

         if (std::abs(delta) < EPS) continue;


         int uu = weight_model.getLaplacianIndex(u, u),
             vv = weight_model.getLaplacianIndex(v, v),
             uv = weight_model.getLaplacianIndex(u, v),
             vu = weight_model.getLaplacianIndex(v, u);

         if (uu == -1 || vv == -1 || uv == -1 || vu == -1) {
             throw std::runtime_error("updateAllEdges: missing matrix entry for edge update");
         }



         weight_model.setEdgeWeight(u, v, new_w);
         weight_model.setEdgeWeight(v, u, new_w);

         // --- Update CSR Laplacian entries ---
         // Diagonal contributions
         m_values[uu] += delta;
         m_values[vv] += delta;

         // Off-diagonal contributions
         m_values[uv] -= delta;
         m_values[vu] -= delta;

     }
 }


void LaplacianSolver::buildLaplacian() {
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


    if (use_dirichlet) {
        m_values_dirichlet = m_values;
        applyDirichletInPlace(m_values_dirichlet);
        hierarchy = std::make_unique<AMG>(std::tie(n, m_row_ptr, m_col_ind, m_values_dirichlet), params);
    } else {
        hierarchy = std::make_unique<AMG>(std::tie(n, m_row_ptr, m_col_ind, m_values), params);
    }
}



std::vector<double> LaplacianSolver::solve(const std::vector<double> &b, double eps) {
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


Eigen::VectorXd LaplacianSolver::solve(const Eigen::VectorXd& b, double eps) {
    const int n_ = b.size();
    std::vector<double> bvec(n_);
    std::memcpy(bvec.data(), b.data(), n_ * sizeof(double));

    result.clear();
    result = this->solve(bvec, eps);
    Eigen::VectorXd eigen_output(result.size());
    std::memcpy(eigen_output.data(), result.data(), n_ * sizeof(double));
    return eigen_output;
}


void LaplacianSolver::updateSolver() {
    // Update hierarchy numeric values (same structure)
    if (use_dirichlet) {
        m_values_dirichlet = m_values;
        applyDirichletInPlace(m_values_dirichlet);
        hierarchy->precond().rebuild(std::tie(n, m_row_ptr, m_col_ind, m_values_dirichlet));
    } else {
        hierarchy->precond().rebuild(std::tie(n, m_row_ptr, m_col_ind, m_values));
    }
}

void LaplacianSolver::applyDirichletInPlace(std::vector<double>& vals) {
     const int r = dirichlet_root;

     // Row r: make it [0 ... 0 1 0 ... 0]
     for (int jj = m_row_ptr[r]; jj < m_row_ptr[r+1]; ++jj) {
         vals[jj] = (m_col_ind[jj] == r) ? 1.0 : 0.0;
     }

     // Column r: zero out A[i,r] for i != r (keep symmetry)
     for (int i = 0; i < n; ++i) {
         if (i == r) continue;
         for (int jj = m_row_ptr[i]; jj < m_row_ptr[i+1]; ++jj) {
             if (m_col_ind[jj] == r) {
                 vals[jj] = 0.0;
                 break;
             }
         }
     }
 }


void LaplacianSolver::setSolverParams(const boost::property_tree::ptree& new_params) {
    params = new_params;
}

void LaplacianSolver::print_params(const boost::property_tree::ptree& prm) {
    boost::property_tree::write_json(std::cout, prm);
}