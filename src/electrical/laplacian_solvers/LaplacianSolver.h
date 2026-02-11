//
// Created by Mert Biyikli on 29.09.25.
//

#ifndef OBLIVIOUSROUTING_LAPLACIANSOLVER_H
#define OBLIVIOUSROUTING_LAPLACIANSOLVER_H

#include "../../datastructures/GraphADJ.h"
#include "../../utils/hash.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <Eigen/Sparse>
#include <Eigen/Dense>

class LaplacianSolver {
public:
    virtual ~LaplacianSolver() = default;

    // virtual void init(std::unordered_map<std::pair<int,int>, double>& edge_weights, int n, bool debug = false) = 0;
    // virtual void buildLaplacian() = 0;

    virtual std::vector<double> solve(const std::vector<double>& b, double eps) = 0;
    virtual Eigen::VectorXd solve(const Eigen::VectorXd& b, double eps) = 0;
    virtual void updateSolver() = 0;
    virtual void buildLaplacian() = 0;
    //virtual void buildLaplacian_(const GraphADJ& g) = 0;
/*
    virtual bool updateEdge(int u, int v, double new_weight) = 0;

    virtual bool checkMatrix() const = 0;
*/
    bool debug = false;
    int n;
    // CSR storage
    std::vector<int> m_row_ptr, m_row_ptr_red;
    std::vector<int> m_col_ind, m_col_ind_red;
    std::vector<double> m_values, m_values_red;

    std::vector<double> result;
    std::vector<double> bvec_buffer, x_buffer;


    std::unordered_map<std::pair<int, int>, double, PairHash> m_edge_weights; // Edge weights

    std::vector<std::vector<double>> adj_edge_weights; // Edge weights

    // handling the indices
    struct indexKey {
        int r, c;
        bool operator==(const indexKey &o) const { return r==o.r && c==o.c; }
    };

    struct indexKeyHash {
        std::size_t operator()(const indexKey &k) const {
            return std::hash<int>()(k.r*73856093 ^ k.c*19349663);
        }
    };
    std::unordered_map<indexKey, int, indexKeyHash> m_indexMap;
    std::unordered_map<std::pair<int,int>, int, PairHash> m_edgeIndexMap;
    int findIndex(int r, int c) const {
        for (int k = m_row_ptr[r]; k < m_row_ptr[r + 1]; ++k)
            if (m_col_ind[k] == c) return k;
        return -1;
    }

    void init(std::vector<double>& _adj_edge_weights,
        int n,
        std::vector<std::pair<int, int>>& edges,
        bool debug = false) {
        this->debug = debug;
        this->n = n;
        m_row_ptr.clear();
        m_col_ind.clear();
        m_values.clear();

        // convert edge list to edge weight map
        m_edge_weights.clear();
        for(int i = 0; i<edges.size(); i++) {
            auto [e1, e2] = edges[i];
            double w = _adj_edge_weights[i];
            m_edge_weights[{e1, e2}] = w;
            m_edge_weights[{e2, e1}] = w; // keep it symmetric
        }

        buildLaplacian();
    }
/*
    void init(std::vector<std::vector<double>>& _adj_edge_weights, int n, const GraphADJ& g, bool debug = false) {
        this->debug = debug;
        this->n = n;
        m_row_ptr.clear();
        m_col_ind.clear();
        m_values.clear();

        // Convert adjacency list to edge weight map
        adj_edge_weights = _adj_edge_weights;

        buildLaplacian_(g);
    }

    void init(std::unordered_map<std::pair<int, int>, double, PairHash>& _edge_weights, int n, bool debug = false) {
        this->debug = debug;
        this->n = n;
        m_row_ptr.clear();
        m_col_ind.clear();
        m_values.clear();

        m_edge_weights = _edge_weights;
        buildLaplacian();
    }*/


    virtual void updateAllEdges(const std::vector<double>& new_weights,
                    const std::vector<std::pair<int,int>>& edges) = 0;
/*
    virtual void updateAllEdges(const std::vector<std::vector<double>>& new_weights, const GraphADJ& G) {
        for (int u = 0; u < n; ++u) {
            for (int idx = 0; idx < G.neighbors(u).size(); ++idx) {
                int v = G.neighbors(u)[idx];

                double weight = new_weights[u][idx];
                double delta = weight - adj_edge_weights[u][idx];
                m_edge_weights[{u, v}] = weight; // Update the weight

                // Update diagonal
                m_values[m_indexMap[{u,u}]] += delta;
                m_values[m_indexMap[{v,v}]] += delta;

                // Update off-diagonals
                m_values[m_indexMap[{u,v}]] -= delta;
                m_values[m_indexMap[{v,u}]] -= delta;

            }
        }
    }

    void updateAllEdges(const std::unordered_map<std::pair<int, int>, double, PairHash>& new_weights) {
        for (const auto& [edge, weight] : new_weights) {

            int u = edge.first;
            int v = edge.second;

            double delta = weight - m_edge_weights[edge];
            m_edge_weights[edge] = weight; // Update the weight

            // Update diagonal
            m_values[m_indexMap[{u,u}]] += delta;
            m_values[m_indexMap[{v,v}]] += delta;

            // Update off-diagonals
            m_values[m_indexMap[{u,v}]] -= delta;
            m_values[m_indexMap[{v,u}]] -= delta;
        }
    }*/

    bool updateEdge(int u, int v, double new_weight) {
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
    bool checkMatrix() const {
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
};

#endif //OBLIVIOUSROUTING_LAPLACIANSOLVER_H