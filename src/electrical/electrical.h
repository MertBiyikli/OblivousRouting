//
// Created by Mert Biyikli on 30.09.25.
//

#ifndef OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H
#define OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H

#include <algorithm>
#include "../solver/solver.h"
#include "../solver/mwu_framework.h"
#include "laplacian_solvers/AMGSolver.h"
#include "../datastructures/IGraph.h"
#include <Eigen/Sparse>


/*
* This is the implementation of the electrical flow based MWU oblivious routing algorithm presented by Goranci et. al. in 2023.
* The idea is to repeatedly invoke an electrical flow computation (Laplacian solve)
* and then update the edge resistances based on the load on the edges. The process is repeated until convergence.
*/
class ElectricalMWU : public LinearObliviousSolverBase, public MWUFramework {
protected:
    // AMG solver instance
    std::unique_ptr<LaplacianSolver> amg;
    int n, m;

    double epsilon = 0.5; // sketching parameter
    double epsilon_L = EPS; // Laplacian solving error
    bool K_initialized = false;
    double K = 1.0; // this is used as an approximation error for the Laplacian Solver( see. Paper for details)

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
    Eigen::SparseMatrix<double> X;          // n × ℓ precomputed RHS
    Eigen::MatrixXd SketchMatrix, SketchMatrix_t; // sketch matrix transposed

public:

    ElectricalMWU(IGraph& g, int root, bool debug = false)
    : LinearObliviousSolverBase(g, root), n(g.getNumNodes()), m(g.getNumEdges()/2) {
        this->debug = debug;
    }

    // entry point
    void computeBasisFlows(LinearRoutingTable &table) override {
        auto start_run = timeNow();
        init(debug);
        run(table);
        scaleFlowDown(table);
        this->solve_time = duration(timeNow() - start_run) - transformation_time;
    }

    virtual void run(LinearRoutingTable &table);
    virtual void scaleFlowDown(LinearRoutingTable &table);
    virtual void getApproxLoad(std::vector<double>& load);
    virtual void init(bool debug = false, boost::property_tree::ptree _params = boost::property_tree::ptree() );
    virtual void initAMGSolver(boost::property_tree::ptree _params);
    virtual void updateEdgeDistances(const std::vector<double>& load);



    // Helpers
    void extractEdges();
    void initEdgeDistances();
    Eigen::SparseMatrix<double> buildIncidence();
    Eigen::MatrixXd getSketchMatrix(int m, int n, double epsilon = 0.5);
    void addFlowToTable(const int& u, Eigen::VectorXd& potential, LinearRoutingTable &table);
};

#endif //OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H