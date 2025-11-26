//
// Created by Mert Biyikli on 22.10.25.
//

#ifndef OBLIVIOUSROUTING_ELECTRICAL_FLOW_PARALLEL_ON_THE_FLY_H
#define OBLIVIOUSROUTING_ELECTRICAL_FLOW_PARALLEL_ON_THE_FLY_H

#include <algorithm>
#include "../datastructures/graph.h"
#include "../utils/hash.h"
#include "../solver/solver.h"
#include "laplacian_solvers/AMGSolverParallel.h"
#include <Eigen/Sparse>
#include "laplacian_solvers/LaplacianSolverFactory.h"

class ElectricalFlowParallelOnTheFly : public ObliviousRoutingSolver{
    private:

    int n, m;
    // LaplacianSolver solver;
    Graph m_graph;
    std::unique_ptr<LaplacianSolver> amg;
    double roh = 0.0;
    double alpha_local = 0.0;
    double inv_m = 0.0;
    std::vector<double> div_accum; // size n
    std::vector<int> tree_parent;
    std::vector<double> edge_weights;             // w_e
    std::vector<double> edge_capacities;          // c_e
    std::vector<double> edge_distances;           // x_e
    std::vector<double> edge_probabilities;       // p_e

    // --- Per-adjacency (directed, aligned with Graph::neighbors(u)) ---
    std::vector<std::vector<int>> adj_edge_ids;   // maps (u, idx) -> edge_id (same id on both directions)


    //std::vector<std::vector<double>> x_edge_distance, p_edge_probability, w_edge_weight, c_edge_capacity;
    double cap_X = 0.0;

    // just for debugging purposes
    std::unordered_map<std::pair<int, int>, double> x_edge2distance, p_edge2probability, w_edges2weights, c_edges2capacities;


    // Preallocate as class members to avoid reallocs
    Eigen::SparseMatrix<double> B; // incidence matrix
    Eigen::SparseMatrix<double> U; // capacity diagonal matrix
    Eigen::SparseMatrix<double> W; // weight diagonal matrix
    Eigen::SparseMatrix<double> SketchMatrix_T; // sketch matrix transposed
    Eigen::SparseMatrix<double> X;          // n × ℓ precomputed RHS
    Eigen::SparseMatrix<double> bw_product_t;

    // Helpers
    void extract_edge_list();
    void initEdgeDistances();
    void updateEdgeDistances(const std::vector<double>& load);
    void buildWeightDiag();
    void refreshWeightMatrix();
    void buildIncidence();

    Eigen::SparseMatrix<double> getSketchMatrix(int m, int n, double epsilon = 0.5);

    public:

    void runSolve(const IGraph& g_) override {
        auto g = dynamic_cast<const Graph&>(g_); // cast to Graph
        this->init(g, debug);
        this->run();
    }
    void storeFlow() override;
    double getFlowForCommodity(int edge_id, int source, int target);

    // std::unordered_map<std::pair<int, int>, std::unordered_map<std::pair<int, int>, double>> f_e_st;
    std::unordered_map<std::pair<int, int>, double> f_e_u; // store the flow of the edge u→x for a fixed vertex x
    int x_fixed;

    // maps (u,v) → edge_id, stores both orientations
    inline long long make_key(int u, int v) const {
        if (u > v) std::swap(u,v); // always store with u<v
        return ( (static_cast<long long>(u) << 32) | static_cast<unsigned int>(v) );
    }

    std::unordered_map<long long, int> edge_id_map;

    // std::vector<WeightedEdge> edges_with_weights; // u<v only


    std::vector<std::pair<int, int> > edges; // u<v only
    // std::vector<double>              weights; // same order
    void run();
    void scaleFlowDown();

    std::vector<double> edge_diffs;
    void getApproxLoad(std::vector<double>& load);
    // Compute the median absolute difference between two node potentials
    // across all ℓ sampled solutions.
    //   diffs_u = potentials[u][*]  (length ℓ)
    //   diffs_v = potentials[v][*]  (length ℓ)
    double recoverNorm(const std::vector<double>& diffs_u,
                                            const std::vector<double>& diffs_v);

    void init(const Graph& g, bool debug = false, const std::string& solver_name = ("amg_parallel"));

/*
    std::vector<int> buildBFSTree(int root);
    void enforceConservation();*/


};

#endif //OBLIVIOUSROUTING_ELECTRICAL_FLOW_PARALLEL_ON_THE_FLY_H