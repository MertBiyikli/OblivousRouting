//
// Created by Mert Biyikli on 09.07.25.
//

#ifndef OBLIVIOUSROUTING_ELECTRICAL_FLOW_NAIVE_H
#define OBLIVIOUSROUTING_ELECTRICAL_FLOW_NAIVE_H

#include <algorithm>
#include "../datastructures/graph.h"
#include "../utils/hash.h"
#include "../solver/solver.h"
#include "laplacian_solvers/AMGSolver.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "../tree_based/frt/raecke_frt_transform.h"
#include "laplacian_solvers/LaplacianSolverFactory.h"

struct FlowPath {
    std::vector<std::pair<int, int>> path;
    double amount;
};


// IMPORTANT TODO: note that the edges are stored only as ordered paris, hence we are essentially missing the rest
//              anti edges...

class ElectricalFlowNaive : public ObliviousRoutingSolver{



    // Preallocate scratch buffers as class members to avoid reallocs
    Eigen::VectorXd scratch_s;
    std::vector<double> scratch_abs;



    std::vector<Eigen::SparseMatrix<double>> Ms;

    Eigen::SparseMatrix<double> W; // weight diagonal matrix
    Eigen::SparseMatrix<double> M_avg;




    int n, m;
    Eigen::SparseMatrix<double> B;
    Eigen::SparseMatrix<double> U;

    Eigen::SparseMatrix<double> bw_product_t;

    double roh = 0.0;
    double alpha_local = 0.0;
    double inv_m = 0.0;
    std::unordered_map<std::pair<int, int>, double, PairHash> x_edge2distance, p_edge2probability, w_edges2weights, c_edges2capacities;
    double cap_X = 0.0;
    //int number_of_iterations = 0;

    // Once‐computed edge list & index mapping

    // TODO: delete this and see if this even influences something...
    std::vector<double>              weights; // same order

    // Helpers
    void extract_edge_list();

    void initEdgeDistances();
    void updateEdgeDistances();
    // std::vector<Eigen::SparseVector<double>>  getRoutingMatrix();
    Eigen::SparseMatrix<double>  getRoutingMatrix();
    void buildWeightDiag();
    void refreshWeightMatrix();


    void buildIncidence();
    Eigen::MatrixXd buildPseudoInverse(Eigen::MatrixXd& L);

    Eigen::SparseMatrix<double> getSketchMatrix(int m, int n, double epsilon = 0.5);
    double recoverNorm(const Eigen::MatrixXd& M, const Eigen::VectorXd& vec);



public:
    // LaplacianSolver solver;
    Graph m_graph;
    std::unique_ptr<LaplacianSolver> amg;

    // std::unordered_map<std::pair<int, int>, std::unordered_map<std::pair<int, int>, double>> f_e_st;
    std::unordered_map<std::pair<int, int>, double, PairHash> f_e_u; // store the flow of the edge u→x for a fixed vertex x
    int x_fixed;
    // store the flow as an weighted adjacency list
    // the first entry denotes the edge and the second
    // entry denotes the flow for the commodity u->x, since x is fixed, we only keep it for u
    // std::vector<std::vector<int>> f_adj_edges;
    // std::vector<std::vector<double>> f_adj_flow;
    std::vector<std::pair<int,int>> edges; // u<v only
    void run();

    std::unordered_map<std::pair<int, int>, double, PairHash> getApproxLoad();
    void init(const Graph& g,bool debug = false,const std::string& solver_name = ("amg_cg"));

    Eigen::SparseMatrix<double> getFinalRoutingMatrix();
    std::unordered_map<std::pair<int, int>, Eigen::VectorXd, PairHash> getRoutingForCommodity(const std::vector<std::pair<int, int> >& _commodity);


    // TODO: remove these
    std::vector<FlowPath> decomposeFlowToPaths(
            int source,
            int sink,
            const Eigen::VectorXd& edge_flow
    );

    std::unordered_map<std::pair<int, int>, std::vector<FlowPath>, PairHash>
    getRoutingPathsForCommodity(const std::vector<std::pair<int, int>>& _commodity);
    double getMaximumCongestion();
    double getCongestion(DemandMap& _demands);
    double _getCongestion(const DemandMap& _demands);

    double getFlowForCommodity(int edge_id, int source, int target);
    double getCongForCommodity(int edge_id, int source, int target);

    void runSolve(const IGraph& g_) override {
        auto g = dynamic_cast<const Graph&>(g_); // cast to Graph
        this->init(g, debug);
        this->run();
    }
    void storeFlow() override;

};

#endif //OBLIVIOUSROUTING_ELECTRICAL_FLOW_NAIVE_H
