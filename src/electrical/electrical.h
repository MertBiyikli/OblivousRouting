//
// Created by Mert Biyikli on 30.09.25.
//

#ifndef OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H
#define OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H

#include <algorithm>
#include "../solver/solver.h"
#include "../solver/mwu_framework.h"
#include "laplacian_solvers/AMGSolver.h"
#include <Eigen/Sparse>

#include "../datastructures/IGraph.h"


/*
* This is the implementation of the algorithm from " Electrical Flows for Polylogarithmic Competitive Oblivious Routing" by Goranci et al.
 */
class ElectricalMWU : public LinearObliviousSolverBase, public MWUFramework {

    // AMG solver instance
    std::unique_ptr<LaplacianSolver> amg;
public:

    ElectricalMWU(IGraph& g, int root, bool debug = false)
    : LinearObliviousSolverBase(g, root), n(g.getNumNodes()), m(g.getNumEdges()/2) {
        this->debug = debug;
    }

    // entry point
    void computeBasisFlows(LinearRoutingTable &table) override {
        init(debug);
        run(table);
        scaleFlowDown(table);
    }

    virtual void run(LinearRoutingTable &table);
    void addFlowToTable(const int& u, Eigen::VectorXd& potential, LinearRoutingTable &table);
    virtual void scaleFlowDown(LinearRoutingTable &table);
    virtual void getApproxLoad(std::vector<double>& load);
    virtual void init(bool debug = false, boost::property_tree::ptree _params = boost::property_tree::ptree() );


    int n, m;
    double K = 1.0; // this is used as an approximation errror for the Laplacian Solver( see. Paper for details)

    double epsilon = 0.5; // sketching parameter
    double epsilon_L = EPS; // Laplacian solving error

    double roh = 0.0;
    double alpha_local = 0.0;
    double inv_m = 0.0;

    double cap_X = 0.0;
    bool debug = false;
    int x_fixed = 0; // fixed node


    std::vector<std::pair<int, int> > edges; // u<v only
    std::vector<double> edge_weights;             // w_e
    std::vector<double> edge_capacities;          // c_e
    std::vector<double> edge_distances;           // x_e
    std::vector<double> edge_probabilities;       // p_e
    std::vector<double> edge_diffs;




    // Preallocate as class members to avoid reallocs
    Eigen::SparseMatrix<double> B; // incidence matrix
    Eigen::SparseMatrix<double> X;          // n × ℓ precomputed RHS

    Eigen::MatrixXd SketchMatrix, SketchMatrix_t; // sketch matrix transposed

    // Helpers
    void extractEdges();
    void initEdgeDistances();
    void buildIncidence();
    virtual void updateEdgeDistances(const std::vector<double>& load);
    virtual Eigen::MatrixXd getSketchMatrix(int m, int n, double epsilon = 0.5);


};

#endif //OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H