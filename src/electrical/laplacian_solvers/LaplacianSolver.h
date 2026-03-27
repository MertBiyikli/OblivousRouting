//
// Created by Mert Biyikli on 29.09.25.
//

#ifndef OBLIVIOUSROUTING_LAPLACIANSOLVER_H
#define OBLIVIOUSROUTING_LAPLACIANSOLVER_H

#include <iostream>

#include "../../datastructures/IGraph.h"
#include "GraphToLaplacian.h"
#include "../../utils/my_math.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

class LaplacianSolver {
protected:
    GraphToLaplacian weight_model;
    bool debug = false;
    int n;
    // CSR storage
    std::vector<int> m_row_ptr;
    std::vector<int> m_col_ind;
    std::vector<double> m_values;

    // reusable buffers
    std::vector<double> result;
    std::vector<double> bvec_buffer, x_buffer;

    bool use_dirichlet = true;
    int dirichlet_root = 0;
    std::vector<double> m_values_dirichlet;

    // this for configuration of the solver, e.g. coarsening and relaxation types
    boost::property_tree::ptree params;

public:
    virtual ~LaplacianSolver() = default;

    virtual std::vector<double> solve(const std::vector<double>& b, double eps) = 0;
    virtual Eigen::VectorXd solve(const Eigen::VectorXd& b, double eps) = 0;
    virtual void updateSolver() = 0;

    virtual void buildLaplacian() = 0;

    void setSolverParams(const boost::property_tree::ptree& new_params) {
        params = new_params;
    }

    void print_params(const boost::property_tree::ptree& prm) {
        boost::property_tree::write_json(std::cout, prm);
    }


    void init(IGraph& g, std::vector<double>& _adj_edge_weights, int n, const std::vector<std::pair<int, int>>& edges, bool debug = false) {
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

    void updateAllEdges(const std::vector<double> &new_weights, const std::vector<std::pair<int, int> > &edges) {
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


    void applyDirichletInPlace(std::vector<double>& vals) {
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
                //ok = false;
            } else if (std::abs(m_values[diag_index]) < 1e-14) {
                std::cerr << "Zero diagonal entry at (" << i << "," << i << ")\n";
                //ok = false;
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

    int findIndex(int r, int c) const {
        for (int k = m_row_ptr[r]; k < m_row_ptr[r + 1]; ++k)
            if (m_col_ind[k] == c) return k;
        return -1;
    }
};

#endif //OBLIVIOUSROUTING_LAPLACIANSOLVER_H