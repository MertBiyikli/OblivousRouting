//
// Created by Mert Biyikli on 09.07.25.
//

#ifndef OBLIVIOUSROUTING_ELECTRICAL_FLOW_NAIVE_H
#define OBLIVIOUSROUTING_ELECTRICAL_FLOW_NAIVE_H

#include <algorithm>
#include "../graph.h"
#include "../utils/hash.h"
#include "laplcian_solver.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>


// IMPORTANT TODO: note that the edges are stored only as ordered paris, hence we are essentially missing the rest
// anti edges...

class ElectricalFlowNaive{
    RaeckeGraph m_graph;
    double roh = 0.0;
    double alpha_local = 0.0;
    double inv_m = 0.0;
    std::unordered_map<std::pair<int, int>, double> x_edge2distance, p_edge2probability, w_edges2weights;
    double cap_X = 0.0;
    int number_of_iterations = 0;

    // Once‐computed edge list & index mapping
    std::vector<std::pair<int,int>> edges; // u<v only
    std::vector<double>              weights; // same order

    // Helpers
    void extract_edge_list();
    void computeM_via_solves(Eigen::SparseMatrix<double> &M_out) {return;}
    void updateEdgeDistances();
    Eigen::SparseMatrix<double> getRoutingMatrix();
    Eigen::SparseMatrix<double> buildWeightDiag(const std::vector<double>&);



    Eigen::SparseMatrix<double> buildIncidence(int V,
                                               std::vector<std::pair<int,int>>& edges);


    Eigen::MatrixXd getSketchMatrix(int m, int n, double epsilon = 0.5);
    double recoverNorm(const Eigen::MatrixXd& M, const Eigen::VectorXd& vec);

public:
    void run();

    std::unordered_map<std::pair<int, int>, double> getApproxLoad();
    void init(const RaeckeGraph& g);


};

#endif //OBLIVIOUSROUTING_ELECTRICAL_FLOW_NAIVE_H
