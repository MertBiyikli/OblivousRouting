//
// Created by Mert Biyikli on 09.07.25.
//

#ifndef OBLIVIOUSROUTING_ELECTRICAL_FLOW_NAIVE_H
#define OBLIVIOUSROUTING_ELECTRICAL_FLOW_NAIVE_H

#include <algorithm>
#include "../graph.h"
#include "../utils/hash.h"
#include "AMGSolver.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "../tree_based/raecke_transform.h"

struct FlowPath {
    std::vector<std::pair<int, int>> path;
    double amount;
};


// IMPORTANT TODO: note that the edges are stored only as ordered paris, hence we are essentially missing the rest
//              anti edges...

class ElectricalFlowNaive{

    std::unordered_map<std::tuple<int, int, int>, double> f_e_st;

    bool debug = false;
    std::vector<Eigen::SparseMatrix<double>> Ms;

    Eigen::SparseMatrix<double> W; // weight diagonal matrix
    Eigen::SparseMatrix<double> M_avg;


    // LaplacianSolver solver;
    Graph m_graph;
    AMGSolver amg;

    int n, m;
    Eigen::SparseMatrix<double> B;
    Eigen::SparseMatrix<double> U;
    double roh = 0.0;
    double alpha_local = 0.0;
    double inv_m = 0.0;
    std::unordered_map<std::pair<int, int>, double> x_edge2distance, p_edge2probability, w_edges2weights, c_edges2capacities;
    double cap_X = 0.0;
    int number_of_iterations = 0;

    // Once‐computed edge list & index mapping

    // TODO: delete this and see if this even influences something...
    std::vector<double>              weights; // same order

    // Helpers
    void extract_edge_list();

    void initEdgeDistances();
    void updateEdgeDistances();
    Eigen::SparseMatrix<double>  getRoutingMatrix();
    void buildWeightDiag();
    void refreshWeightMatrix();


    void buildIncidence();
    Eigen::MatrixXd buildPseudoInverse(Eigen::MatrixXd& L);

    Eigen::SparseMatrix<double> getSketchMatrix(int m, int n, double epsilon = 0.5);
    double recoverNorm(const Eigen::MatrixXd& M, const Eigen::VectorXd& vec);

public:
    std::unordered_map<std::pair<int, int>, double> f_e_u; // store the flow of the edge u→x for a fixed vertex x
    int x_fixed;
    // store the flow as an weighted adjacency list
    // the first entry denotes the edge and the second
    // entry denotes the flow for the commodity u->x, since x is fixed, we only keep it for u
    // std::vector<std::vector<int>> f_adj_edges;
    // std::vector<std::vector<double>> f_adj_flow;
    std::vector<std::pair<int,int>> edges; // u<v only
    void run();

    std::unordered_map<std::pair<int, int>, double> getApproxLoad();
    void init(const Graph& g, bool debug = false);

    Eigen::SparseMatrix<double> getFinalRoutingMatrix();
    std::unordered_map<std::pair<int, int>, Eigen::VectorXd> getRoutingForCommodity(const std::vector<std::pair<int, int> >& _commodity);


    // TODO: remove these
    std::vector<FlowPath> decomposeFlowToPaths(
            int source,
            int sink,
            const Eigen::VectorXd& edge_flow
    );

    std::unordered_map<std::pair<int, int>, std::vector<FlowPath>>
    getRoutingPathsForCommodity(const std::vector<std::pair<int, int>>& _commodity);
    double getMaximumCongestion();
    double getCongestion(DemandMap& _demands);
};

#endif //OBLIVIOUSROUTING_ELECTRICAL_FLOW_NAIVE_H
