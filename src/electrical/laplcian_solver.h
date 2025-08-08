//
// Created by Mert Biyikli on 09.07.25.
//

#ifndef OBLIVIOUSROUTING_LAPLCIAN_SOLVER_H
#define OBLIVIOUSROUTING_LAPLCIAN_SOLVER_H


#include <julia.h>
#include <julia_version.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <julia/julia.h>
#include "../graph.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unordered_map>
#include "../utils/hash.h"

JULIA_DEFINE_FAST_TLS // only define this once, in an executable (not in a shared library) if you want fast code.


/* TODO:
 *  1 Major Todos to boost running time:
 *   - dynamically adjust the weights of the adjacency matrix instead of creating a new one each time
 */
class LaplacianSolver {
    bool debug = false; // Debug flag for printing debug information
    RaeckeGraph m_graph;
    std::unordered_map<std::pair<int, int>, double> w_edges2weights; // Edge weights

    bool initted;
    int n, m;

    // storing the adjacency matrix
    std::vector<int32_t> colptr, rowval;
    std::vector<double> nzval;

    Eigen::SparseMatrix<double> L;

    // store the demand vector
    Eigen::VectorXd bvec; // Demand vector for the linear system

    // store julia vector to avoid recomputing everything from scratch
    jl_array_t* colptr_julia;
    jl_array_t* rowval_julia;
    jl_array_t* nzval_julia;
    jl_array_t* bvec_julia; // Demand vector for the linear system in Julia

    jl_value_t* arg1;
    jl_value_t* arg2;

    // initialize the Julia arrays and setting the adjacency matrix
    void buildCSC();
public:
    LaplacianSolver();
    Eigen::VectorXd solve(const Eigen::VectorXd &b);

    void init(const RaeckeGraph &g, const std::unordered_map<std::pair<int, int>, double> &_w_edges2weights);

    void updateEdgeWeights(std::pair<int, int> edge, double weight);

    void NonSyncUpdateEdgeWeights(std::pair<int, int> edge, double weight);
    void SyncEdgeWeights();
    void setDebug(bool d) {
        debug = d;
    }

    Eigen::SparseMatrix<double> GetLaplacianMatrix() const {
        return L;
    }
};

#endif //OBLIVIOUSROUTING_LAPLCIAN_SOLVER_H
