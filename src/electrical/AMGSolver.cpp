//
// Created by Mert Biyikli on 06.08.25.
//

#include "AMGSolver.h"

void AMGSolver::init(std::unordered_map<std::pair<int, int>, double>& _edge_weights, int n, bool debug) {
    this->debug = debug;
    this->n = n;
    m_row_ptr.clear();
    m_col_ind.clear();
    m_values.clear();

    m_edge_weights = _edge_weights;
    buildLaplacian();
}

void AMGSolver::buildLaplacian() {
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
}



std::vector<double> AMGSolver::solve(const std::vector<double> &b) {
    if (b.size() != n)
        throw std::runtime_error("RHS size mismatch");

    auto bvec= b;
    double mean_bvec = std::accumulate(bvec.begin(), bvec.end(), 0.0) / n;
    // divide each entry of bvec by the mean to ensure zero-sum
    for(auto& val : bvec) {
        val -= mean_bvec;
    }

    Solver::params prm;
    prm.solver.tol = 1e-12;

    prm.precond.direct_coarse = false;

    if(debug) {
        std::cout << "row_ptr\n";
        for (auto r: m_row_ptr) std::cout << r << " ";
        std::cout << "\nm_col_ind\n";
        for (auto c: m_col_ind) std::cout << c << " ";
        std::cout << "\nm_values\n";
        for (auto v: m_values) std::cout << v << " ";
        std::cout << std::endl;
    }


    // Construct solver
    Solver solver(std::tie(n, m_row_ptr, m_col_ind, m_values), prm);

    std::vector<double> x(n, 0.0);
    auto res = solver(bvec, x);
    if(debug) {
        std::cout << "[AMGSolver] Iterations: " << std::get<0>(res)
                  << "  Error: " << std::get<1>(res) << std::endl;
    }
    // --- Enforce zero-mean solution (orthogonal to nullspace) ---
    double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / n;
    for(auto& val : x) {
        val -= mean_x; // Project solution to zero-mean
    }

    return x;
}

Eigen::VectorXd AMGSolver::solve(const Eigen::VectorXd& b) {
    std::vector<double> bvec;
    for(int i = 0; i < b.size(); ++i) {
        bvec.emplace_back(b[i]);
    }

    auto res = this->solve(bvec);
    Eigen::VectorXd result(res.size());
    for (int i = 0; i < res.size(); ++i) {
        result[i] = res[i];
    }
    return result;
}


bool AMGSolver::updateEdge(int u, int v, double new_weight) {
    auto it = m_edge_weights.find({u, v});
    if (it == m_edge_weights.end()) {
        return false;
    }else{
        double delta = new_weight - it->second;
        it->second = new_weight; // Update the weight

        // Update diagonal
        m_values[m_indexMap[{u,u}]] += delta;
        m_values[m_indexMap[{v,v}]] += delta;

        // Update off-diagonals
        m_values[m_indexMap[{u,v}]] -= delta;
        m_values[m_indexMap[{v,u}]] -= delta;

        return true;
    }
}

// TODO : implement this
void AMGSolver::updateSolver() {
    // Update hierarchy numeric values (same structure)
    // solver.precond().rebuild(std::tie(n, m_row_ptr, m_col_ind, m_values));

}

void AMGSolver::pinVertex(int v) {
    for (int idx = m_row_ptr[v]; idx < m_row_ptr[v+1]; ++idx)
        m_values[idx] = 0.0;
    for (int row = 0; row < n; ++row) {
        for (int idx = m_row_ptr[row]; idx < m_row_ptr[row+1]; ++idx) {
            if (m_col_ind[idx] == v)
                m_values[idx] = 0.0;
        }
    }
    int diag = m_indexMap[{v,v}];
    m_values[diag] = 1.0;
}

bool AMGSolver::checkMatrix() const {
    bool ok = true;

    // --- 1. Check for symmetry and duplicates ---
    for (int row = 0; row < n; ++row) {
        std::unordered_set<int> cols_in_row;
        for (int idx = m_row_ptr[row]; idx < m_row_ptr[row + 1]; ++idx) {
            int col = m_col_ind[idx];
            double val = m_values[idx];

            // Duplicates check:
            if (!cols_in_row.insert(col).second) {
                std::cerr << "Duplicate entry at row " << row << ", col " << col << "\n";
                ok = false;
            }

            // Symmetry check:
            if (col != row) {
                int rev_index = findIndex(col, row);
                if (rev_index == -1) {
                    std::cerr << "Missing symmetric entry at (" << col << "," << row << ")\n";
                    ok = false;
                } else if (std::abs(m_values[rev_index] - val) > 1e-12) {
                    std::cerr << "Matrix not symmetric at (" << row << "," << col
                              << ") vs (" << col << "," << row << ")\n";
                    ok = false;
                }
            }
        }
    }

    // --- 2. Check diagonal entries ---
    for (int i = 0; i < n; ++i) {
        int diag_index = findIndex(i, i);
        if (diag_index == -1) {
            std::cerr << "Missing diagonal entry at (" << i << "," << i << ")\n";
            ok = false;
        } else if (std::abs(m_values[diag_index]) < 1e-14) {
            std::cerr << "Zero diagonal entry at (" << i << "," << i << ")\n";
            ok = false;
        }
    }

    // --- 3. Check row sums (Laplacian property) ---
    for (int row = 0; row < n; ++row) {
        double sum = 0.0;
        for (int idx = m_row_ptr[row]; idx < m_row_ptr[row + 1]; ++idx) {
            sum += m_values[idx];
        }
        if (std::abs(sum) > 1e-10) {
            std::cerr << "Row " << row << " does not sum to 0: sum = " << sum << "\n";
            // not necessarily fatal, but unusual for Laplacian
        }
    }

    return ok;
}
