//
// Created by Mert Biyikli on 29.09.25.
//

#include "AMGBiCGStab.h"
#include <amgcl/adapter/crs_tuple.hpp>

std::vector<double> AMGBiCGStabSolver::solve(const std::vector<double>& b) {
    if (!solver) throw std::runtime_error("Solver not initialized");
    if (b.size() != n) throw std::runtime_error("RHS size mismatch");

    auto bvec = b;
    double mean = std::accumulate(bvec.begin(), bvec.end(), 0.0) / n;
    for (auto& val : bvec) val -= mean;

    std::vector<double> x(n, 0.0);
    auto res = (*solver)(bvec, x);

    if (debug) {
        std::cout << "[AMGBiCGStabSolver] Iter: " << std::get<0>(res)
                  << " Err: " << std::get<1>(res) << std::endl;
    }

    // enforce zero mean
    double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / n;
    for (auto& val : x) val -= mean_x;
    return x;
}

Eigen::VectorXd AMGBiCGStabSolver::solve(const Eigen::VectorXd& b) {
    std::vector<double> bvec(b.data(), b.data() + b.size());
    auto res = solve(bvec);
    Eigen::VectorXd result(res.size());
    for (int i = 0; i < res.size(); ++i) result[i] = res[i];
    return result;
}

void AMGBiCGStabSolver::updateSolver() {
    return;
}

void AMGBiCGStabSolver::buildLaplacian() {
    // --- Step 1: Aggregate edge contributions into a Laplacian map ---
    struct pair_hash {
        size_t operator()(const std::pair<int,int>& p) const {
            return std::hash<int>()(p.first * 73856093 ^ p.second * 19349663);
        }
    };
    std::unordered_map<std::pair<int,int>, double, pair_hash> L;

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

    if(debug){
        if (!checkMatrix()) {
            throw std::runtime_error("Matrix structure is invalid");
        }
    }
    Solver::params prm;
    prm.solver.tol = 1e-12;
    solver = std::make_unique<Solver>(std::tie(n, m_row_ptr, m_col_ind, m_values), prm);
}

void AMGBiCGStabSolver::buildLaplacian_(const Graph &g) {
    return;
}
