//
// Created by Mert Biyikli on 06.08.25.
//

#include "AMGSolver.h"
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

void AMGSolver::check_openmp_runtime() {
#ifdef _OPENMP
    std::cout << "✅ OpenMP is enabled at compile time.\n";
    std::cout << "🧵 omp_get_num_threads(): ";
#pragma omp parallel
    {
#pragma omp master
        {
            std::cout << omp_get_num_threads() << "\n";
        }
    }
#endif
}/*
void AMGSolver::init(std::unordered_map<std::pair<int, int>, double>& _edge_weights, int n, bool debug) {
    this->debug = debug;
    this->n = n;
    m_row_ptr.clear();
    m_col_ind.clear();
    m_values.clear();

    m_edge_weights = _edge_weights;
    buildLaplacian();
}*/
/*
void AMGSolver::buildLaplacian_(const GraphADJ& g) {
    check_openmp_runtime();
    // --- Step 1: Aggregate edge contributions into a Laplacian map ---
    struct pair_hash {
        size_t operator()(const std::pair<int,int>& p) const {
            return std::hash<int>()(p.first * 73856093 ^ p.second * 19349663);
        }
    };
    std::unordered_map<std::pair<int,int>, double, pair_hash> L;

    for (int u = 0; u < n; ++u) {
        for (int idx = 0; idx<g.neighbors(u).size(); idx++) {
            int v = g.neighbors(u)[idx];
            double w = adj_edge_weights[u][idx];
            if (u == v) continue;

            // Laplacian contributions:
            L[{u,u}] += w;
            L[{v,v}] += w;
            L[{u,v}] -= w;
            L[{v,u}] -= w;
        }
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
    // solver = std::make_unique<Solver>(std::tie(n, m_row_ptr, m_col_ind, m_values), prm);
    hierarchy = std::make_unique<AMG>(std::tie(n, m_row_ptr, m_col_ind, m_values));
}*/


void AMGSolver::buildLaplacian() {
    check_openmp_runtime();
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
    if (debug ) {
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


    if(debug) {
        std::cout << "row_ptr\n";
        for (auto r: m_row_ptr) std::cout << r << " ";
        std::cout << "\nm_col_ind\n";
        for (auto c: m_col_ind) std::cout << c << " ";
        std::cout << "\nm_values\n";
        for (auto v: m_values) std::cout << v << " ";
        std::cout << std::endl;
    }

    assert(checkMatrix() && "Matrix structure is invalid during buildLaplacian");

    hierarchy = std::make_unique<AMG>(std::tie(n, m_row_ptr, m_col_ind, m_values));

    // Construct solver
    // solver = std::make_unique<Solver>(std::tie(n, m_row_ptr, m_col_ind, m_values), prm);

}



std::vector<double> AMGSolver::solve(const std::vector<double> &b, double eps) {
    if (b.size() != n)
        throw std::runtime_error("RHS size mismatch");

    x_buffer.assign(n, 0.0);
    bvec_buffer.assign(b.begin(), b.end());

    // --- Enforce zero-mean RHS (orthogonal to nullspace) ---
    double mean_bvec = std::reduce(bvec_buffer.begin(), bvec_buffer.end(), 0.0) / n;
    // divide each entry of bvec by the mean to ensure zero-sum
    for(auto& val : bvec_buffer) {
        val -= mean_bvec;
    }

    amgcl::solver::cg<Backend>::params solver_params;
    solver_params.tol = eps;

    amgcl::solver::cg<Backend> solver(n,solver_params);

    if (true) {
        if(debug) {
            const auto res = solver(*hierarchy, bvec_buffer, x_buffer);
            std::cout << "[AMGSolver] Iterations: " << std::get<0>(res)
                      << "  Error: " << std::get<1>(res) << std::endl;
        }else {
            solver(*hierarchy, bvec_buffer, x_buffer);
        }
        // --- Enforce zero-mean solution (orthogonal to nullspace) ---
        double mean_x = std::reduce( x_buffer.begin(), x_buffer.end()) / n;
        for(auto& val : x_buffer) {
            val -= mean_x; // Project solution to zero-mean
        }
    }else {
        std::cout << "[AMGSolver] No solver found." << std::endl;
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

/*
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
}*/

void AMGSolver::updateSolver() {
/* BUG: In some iterations the AMGSolver outputs "Zero sum in skyline_lu factorization"
 *      The error is putput when the factorize method is called.
 *
*/
    assert(checkMatrix() && "Matrix structure is invalid during updateSolver");
    // Update hierarchy numeric values (same structure)
    hierarchy->rebuild(std::tie(n, m_row_ptr, m_col_ind, m_values));
}

