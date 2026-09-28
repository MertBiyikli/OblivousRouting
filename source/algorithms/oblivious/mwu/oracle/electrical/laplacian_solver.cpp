//
// Created by Mert Biyikli on 20.03.26.
//

#include "algorithms/oblivious/mwu/oracle/electrical/laplacian_solver.h"
#include "utils/my_math.h"
#include "utils/hash.h"
#include <amgcl/backend/builtin.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include "amgcl/solver/runtime.hpp"
#include <amgcl/make_solver.hpp>
#include <cassert>
#include <stdexcept>

void LaplacianSolver::init(optimized::Graph<EdgeData>& g, std::vector<double>& _adj_edge_weights, int n, const std::vector<std::pair<int, int>>& edges, bool debug) {
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

void LaplacianSolver::init(
    int n_,
    const std::vector<std::pair<int, int>>& edges,
    const std::vector<double>& edge_weights,
    bool debug_
) {
    n = n_;
    debug = debug_;

    if (static_cast<int>(edge_weights.size()) != static_cast<int>(edges.size())) {
        throw std::runtime_error("LaplacianSolver::init: edge_weights and edges size mismatch");
    }

    m_values = edge_weights;

    // If GraphToLaplacian requires an IGraph, skip it here and build directly.
    m_row_ptr.assign(n + 1, 0);

    for (const auto& [u, v] : edges) {
        ++m_row_ptr[u + 1];
        ++m_row_ptr[v + 1];
    }

    for (int i = 1; i <= n; ++i) {
        m_row_ptr[i] += m_row_ptr[i - 1];
    }

    m_col_ind.assign(m_row_ptr.back(), 0);
    m_values.assign(m_row_ptr.back(), 0.0);

    std::vector<int> cursor = m_row_ptr;

    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto [u, v] = edges[e];
        const double w = edge_weights[e];

        const int pu = cursor[u]++;
        m_col_ind[pu] = v;
        m_values[pu] = -w;

        const int pv = cursor[v]++;
        m_col_ind[pv] = u;
        m_values[pv] = -w;
    }

    // Add diagonal entries.
    std::vector<double> diag(n, 0.0);
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto [u, v] = edges[e];
        const double w = edge_weights[e];

        diag[u] += w;
        diag[v] += w;
    }

    // Simpler: rebuild CSR with diagonal included.
    std::vector<std::vector<std::pair<int, double>>> rows(n);

    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto [u, v] = edges[e];
        const double w = edge_weights[e];

        rows[u].push_back({v, -w});
        rows[v].push_back({u, -w});
    }

    for (int u = 0; u < n; ++u) {
        rows[u].push_back({u, diag[u]});
    }

    m_row_ptr.assign(n + 1, 0);
    for (int u = 0; u < n; ++u) {
        m_row_ptr[u + 1] = m_row_ptr[u] + static_cast<int>(rows[u].size());
    }

    m_col_ind.clear();
    m_values.clear();

    m_col_ind.reserve(m_row_ptr.back());
    m_values.reserve(m_row_ptr.back());

    for (int u = 0; u < n; ++u) {
        for (const auto& [v, value] : rows[u]) {
            m_col_ind.push_back(v);
            m_values.push_back(value);
        }
    }

    if (use_dirichlet) {
        applyDirichletInPlace(m_values);
    }

    weight_model.init(edges, edge_weights, n);
    for (int e = 0; e < edges.size(); e++) {
        int u = edges[e].first;
        int v = edges[e].second;
        double w = edge_weights[e];
        weight_model.setEdgeWeight(u, v, w);
        weight_model.setEdgeWeight(v, u, w); // keep it symmetric
    }
    buildLaplacian();
    //updateSolver();
}

Result<void> LaplacianSolver::initGrounded(const int node_count,const std::vector<std::pair<int, int>>& edges,const std::vector<double>& edge_weights,const int grounded_vertex) {
    hierarchy.reset();

    if (node_count <= 0) {
        return makeErrorMessage(
            ErrorCode::InvalidGraph,
            "Grounded Laplacian requires at least one vertex."
        );
    }

    if (grounded_vertex < 0 ||
        grounded_vertex >= node_count) {
        return makeErrorMessage(
            ErrorCode::InvalidGraph,
            "Invalid grounded vertex."
        );
    }

    if (edges.size() != edge_weights.size()) {
        return makeErrorMessage(
            ErrorCode::InvalidGraph,
            "Local edges and conductances have different sizes."
        );
    }

    n = node_count;
    use_dirichlet = true;
    dirichlet_root = grounded_vertex;

    /*
     * Validate connectivity. Grounding one vertex only removes
     * the nullspace when the local graph is connected.
     */
    std::vector<std::vector<int>> adjacency(n);

    for (std::size_t edge_id = 0;
         edge_id < edges.size();
         ++edge_id) {
        const auto [u, v] = edges[edge_id];
        const double weight = edge_weights[edge_id];

        if (u < 0 || v < 0 ||
            u >= n || v >= n) {
            return makeErrorMessage(
                ErrorCode::InvalidGraph,
                "Local edge has an invalid endpoint."
            );
        }

        if (u == v) {
            return makeErrorMessage(
                ErrorCode::InvalidGraph,
                "Local electrical graph contains a self-loop."
            );
        }

        if (!std::isfinite(weight) ||
            weight <= 0.0) {
            return makeErrorMessage(
                ErrorCode::InvalidGraph,
                "Local conductance must be finite and positive."
            );
        }

        adjacency[u].push_back(v);
        adjacency[v].push_back(u);
    }

    if (n > 1) {
        std::vector<bool> visited(n, false);
        std::queue<int> queue;

        visited[grounded_vertex] = true;
        queue.push(grounded_vertex);

        int visited_count = 0;

        while (!queue.empty()) {
            const int u = queue.front();
            queue.pop();

            ++visited_count;

            for (const int v : adjacency[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    queue.push(v);
                }
            }
        }

        if (visited_count != n) {
            return makeErrorMessage(
                ErrorCode::InvalidGraph,
                "Local cluster graph is disconnected. Visited " +
                    std::to_string(visited_count) +
                    " of " +
                    std::to_string(n) +
                    " vertices."
            );
        }
    }

    /*
     * Build the local Laplacian directly.
     *
     * Do not use buildLaplacian() here, because that method assumes
     * that both directed orientations are present.
     */
    std::vector<std::unordered_map<int, double>> rows(n);

    for (std::size_t edge_id = 0;
         edge_id < edges.size();
         ++edge_id) {
        const auto [u, v] = edges[edge_id];
        const double weight = edge_weights[edge_id];

        rows[u][u] += weight;
        rows[v][v] += weight;

        rows[u][v] -= weight;
        rows[v][u] -= weight;
    }

    /*
     * Every row must have an explicit diagonal.
     */
    for (int vertex = 0; vertex < n; ++vertex) {
        if (!rows[vertex].contains(vertex)) {
            rows[vertex][vertex] = 0.0;
        }
    }

    m_row_ptr.assign(n + 1, 0);

    for (int row = 0; row < n; ++row) {
        m_row_ptr[row + 1] =
            m_row_ptr[row] +
            static_cast<int>(rows[row].size());
    }

    m_col_ind.clear();
    m_values.clear();

    m_col_ind.reserve(m_row_ptr.back());
    m_values.reserve(m_row_ptr.back());

    for (int row = 0; row < n; ++row) {
        std::vector<std::pair<int, double>> entries(
            rows[row].begin(),
            rows[row].end()
        );

        std::sort(
            entries.begin(),
            entries.end(),
            [](const auto& lhs, const auto& rhs) {
                return lhs.first < rhs.first;
            }
        );

        for (const auto& [column, value] : entries) {
            if (!std::isfinite(value)) {
                return makeErrorMessage(
                    ErrorCode::InvalidGraph,
                    "Local Laplacian contains a non-finite entry."
                );
            }

            m_col_ind.push_back(column);
            m_values.push_back(value);
        }
    }

    /*
     * Impose x[grounded_vertex] = 0:
     *
     * - grounded row becomes the identity row;
     * - grounded column is zeroed;
     * - symmetry and positive definiteness are preserved.
     */
    m_values_dirichlet = m_values;
    applyDirichletInPlace(m_values_dirichlet);
    boost::property_tree::ptree solver_params = params;

    /*
     * The grounded Laplacian is symmetric positive definite, so CG
     * is an appropriate outer solver.
     */
    solver_params.put("solver.type","cg");
    solver_params.put("solver.tol",1e-10);
    solver_params.put("solver.maxiter",1000);
    solver_params.put("precond.coarsening.type","smoothed_aggregation");
    solver_params.put("precond.relax.type","spai0");
    try {
        hierarchy = std::make_unique<AMG>(
            std::tie(
                n,
                m_row_ptr,
                m_col_ind,
                m_values_dirichlet
            ),
            params
        );
    } catch (const std::exception& error) {
        return makeErrorMessage(
            ErrorCode::SolverFailed,
            "AMGCL hierarchy construction failed: " +
                std::string(error.what())
        );
    }

    if (!hierarchy) {
        return makeErrorMessage(
            ErrorCode::SolverFailed,
            "AMGCL hierarchy was not constructed."
        );
    }

    return {};
}

Result<void> LaplacianSolver::updateAllEdges(const std::vector<double> &new_weights, const std::vector<std::pair<int, int> > &edges) {
     assert(new_weights.size() == edges.size());

     for (size_t e = 0; e < edges.size(); ++e) {
         int u = edges[e].first;
         int v = edges[e].second;
         double old_w = weight_model.getEdgeWeight(u, v);
         double new_w = std::max(new_weights[e], EPS);
         double delta = new_w - old_w;

         if (std::abs(delta) < EPS) continue;


         auto uu = weight_model.getLaplacianIndex(u, u),
             vv = weight_model.getLaplacianIndex(v, v),
             uv = weight_model.getLaplacianIndex(u, v),
             vu = weight_model.getLaplacianIndex(v, u);

         if (!uu
             || !vv
             || !uv
             || !vu) {
             return getError(uu);
         }



         weight_model.setEdgeWeight(u, v, new_w);
         weight_model.setEdgeWeight(v, u, new_w);

         // --- Update CSR Laplacian entries ---
         // Diagonal contributions
         m_values[uu.value()] += delta;
         m_values[vv.value()] += delta;

         // Off-diagonal contributions
         m_values[uv.value()] -= delta;
         m_values[vu.value()] -= delta;

     }
    return {};
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



Result<std::vector<double>>
LaplacianSolver::solve(
    const std::vector<double>& b,
    const double /* eps */
) {
    if (b.size() != static_cast<std::size_t>(n)) {
        return makeErrorMessage(
            ErrorCode::InvalidDemand,
            "AMGCL RHS size mismatch."
        );
    }

    if (!hierarchy) {
        return makeErrorMessage(
            ErrorCode::SolverFailed,
            "AMGCL solver is not initialized."
        );
    }

    std::vector<double> rhs = b;
    std::vector<double> potential(n, 0.0);

    if (use_dirichlet) {
        if (dirichlet_root < 0 ||
            dirichlet_root >= n) {
            return makeErrorMessage(
                ErrorCode::InvalidGraph,
                "Invalid Dirichlet root."
            );
        }

        rhs[dirichlet_root] = 0.0;
        potential[dirichlet_root] = 0.0;
    } else {
        const double mean =
            std::reduce(
                rhs.begin(),
                rhs.end(),
                0.0
            ) /
            static_cast<double>(n);

        for (double& value : rhs) {
            value -= mean;
        }
    }

    for (const double value : rhs) {
        if (!std::isfinite(value)) {
            return makeErrorMessage(
                ErrorCode::InvalidDemand,
                "AMGCL RHS contains a non-finite value."
            );
        }
    }

    std::size_t iterations = 0;
    double relative_error = 0.0;

    try {
        /*
         * hierarchy is already:
         *
         *     iterative solver + AMG preconditioner
         *
         * Do not wrap another CG solver around it.
         */
        const auto solve_result =
            (*hierarchy)(rhs, potential);

        iterations =
            std::get<0>(solve_result);

        relative_error =
            std::get<1>(solve_result);
    } catch (const std::exception& error) {
        return makeErrorMessage(
            ErrorCode::SolverFailed,
            "AMGCL solve failed: " +
                std::string(error.what())
        );
    }

    if (debug) {
        std::cout
            << "[AMGCL] n=" << n
            << " iterations=" << iterations
            << " relative_error=" << relative_error
            << '\n';
    }

    if (!std::isfinite(relative_error)) {
        return makeErrorMessage(
            ErrorCode::SolverFailed,
            "AMGCL returned a non-finite residual."
        );
    }

    if (!allFinite(potential)) {
        return makeErrorMessage(
            ErrorCode::SolverFailed,
            "AMGCL returned non-finite potentials. "
            "Iterations=" +
                std::to_string(iterations) +
                ", relative residual=" +
                std::to_string(relative_error)
        );
    }

    if (use_dirichlet) {
        potential[dirichlet_root] = 0.0;
    }

    return potential;
}

bool LaplacianSolver::allFinite(const std::vector<double>& vec) {
    for (const double& v : vec) {
        if (!std::isfinite(v)) {
            return false;
        }
    }
    return true;
}


Result<Eigen::VectorXd> LaplacianSolver::solve(const Eigen::VectorXd& b, double eps) {
    const int n_ = b.size();
    std::vector<double> bvec(n_);
    std::memcpy(bvec.data(), b.data(), n_ * sizeof(double));

    Result<std::vector<double>> result;
    result = this->solve(bvec, eps);
    if (!result) {
        return getError(result);
    }

    Eigen::VectorXd eigen_output(result.value().size());
    std::memcpy(eigen_output.data(), result.value().data(), n_ * sizeof(double));

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