//
// Created by Mert Biyikli on 06.08.25.
//

#ifndef OBLIVIOUSROUTING_AMGSOLVER_H
#define OBLIVIOUSROUTING_AMGSOLVER_H

#include <vector>
#include <tuple>
#include <iostream>
#include <unordered_set>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_builder.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
//#include <amgcl/adapter/crs.hpp>



#include "../graph.h"
#include <unordered_map>
#include "../utils/hash.h"

class AMGSolver{
private:
    bool debug = false;
    // AMG Solver define types
    typedef amgcl::make_solver<
            amgcl::amg<
                    amgcl::backend::builtin<double>,
                    amgcl::coarsening::smoothed_aggregation,
                    amgcl::relaxation::spai0
            >,
            amgcl::solver::cg<amgcl::backend::builtin<double>>
    > Solver;



    int n;

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
    std::unordered_map<std::pair<int,int>, int> m_edgeIndexMap;
    int findIndex(int r, int c) const {
        for (int k = m_row_ptr[r]; k < m_row_ptr[r + 1]; ++k)
            if (m_col_ind[k] == c) return k;
        return -1;
    }

public:

    // CSR storage
    std::vector<int> m_row_ptr, m_row_ptr_red;
    std::vector<int> m_col_ind, m_col_ind_red;
    std::vector<double> m_values, m_values_red;

    int n_reduced = 0;

    std::unordered_map<std::pair<int, int>, double> m_edge_weights; // Edge weights
    void init(std::unordered_map<std::pair<int, int>, double>& _edge_weights, int n, bool debug = false);
    void buildLaplacian();
    // void buildReducedLaplacian();
    std::vector<double> solve(const std::vector<double> &b);
    Eigen::VectorXd solve(const Eigen::VectorXd &b);
    // Eigen::VectorXd solve_reduced(const Eigen::VectorXd &b_full);

    void updateSolver();
    bool updateEdge(int u, int v, double new_weight);

    bool checkMatrix() const;
    void pinVertex(int v);

};

#endif //OBLIVIOUSROUTING_AMGSOLVER_H
