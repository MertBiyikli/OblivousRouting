//
// Created by Mert Biyikli on 30.09.25.
//

#ifndef OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H
#define OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H

#include <algorithm>
#include "../datastructures/GraphADJ.h"
#include "../utils/hash.h"
#include "../solver/solver.h"
#include "../solver/mwu_framework.h"
#include "laplacian_solvers/AMGSolver.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "../datastructures/GraphCSR.h"
// #include "../tree_based/frt/raecke_frt_transform.h"
#include "laplacian_solvers/LaplacianSolverFactory.h"

struct WeightedEdge {
    std::pair<int, int> edge; // (u,v) with u<v
    double weight;
};

class ElectricalFlowOptimized : public LinearObliviousSolverBase, public MWUFramework {
    private:

    int n, m;
    // LaplacianSolver solver;
    // const GraphCSR* m_graph = nullptr;
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


    //std::vector<std::vector<double>> x_edge_distance, p_edge_probability, w_edge_weight, c_edge_capacity;
    double cap_X = 0.0;
    // int iteration_count = 0;
    bool debug = false;
    int x_fixed = 0; // fixed node

    // just for debugging purposes
    std::unordered_map<std::pair<int, int>, double, PairHash> x_edge2distance, p_edge2probability, w_edges2weights, c_edges2capacities;


    // Preallocate as class members to avoid reallocs
    Eigen::SparseMatrix<double> B; // incidence matrix
    Eigen::SparseMatrix<double> U; // capacity diagonal matrix
    Eigen::SparseMatrix<double> W; // weight diagonal matrix
    Eigen::SparseMatrix<double> SketchMatrix_T; // sketch matrix transposed
    Eigen::SparseMatrix<double> X;          // n × ℓ precomputed RHS
    Eigen::SparseMatrix<double> bw_product_t;

    // Helpers
    void extract_edge_list(); // TODO: refactor this
    void initEdgeDistances();
    void updateEdgeDistances(const std::vector<double>& load);
    void buildWeightDiag();
    void refreshWeightMatrix();
    void buildIncidence();

    Eigen::SparseMatrix<double> getSketchMatrix(int m, int n, double epsilon = 0.5);

    public:

    ElectricalFlowOptimized(IGraph& g, int root, bool debug = false):LinearObliviousSolverBase(g, root) {}

    void computeBasisFlows(LinearRoutingTable &table) override {
        init(debug, "amg_cg");

        run(table);
        scaleFlowDown(table);

    }




    std::vector<std::pair<int, int> > edges; // u<v only
    // std::vector<double>              weights; // same order
    void run(LinearRoutingTable &table);
    void scaleFlowDown(LinearRoutingTable &table);

    std::vector<double> edge_diffs;
    void getApproxLoad(std::vector<double>& load);
    // Compute the median absolute difference between two node potentials
    // across all ℓ sampled solutions.
    //   diffs_u = potentials[u][*]  (length ℓ)
    //   diffs_v = potentials[v][*]  (length ℓ)
    double recoverNorm(const std::vector<double>& diffs_u,
                                            const std::vector<double>& diffs_v);

    void init( bool debug = false, const std::string& solver_name = ("amg_cg"));

/*
    std::vector<int> buildBFSTree(int root);
    void enforceConservation();*/


};

#endif //OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H