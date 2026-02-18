//
// Created by Mert Biyikli on 04.02.26.
//

#ifndef OBLIVIOUSROUTING_ELECTRICAL_FLOW_PAR_BATCHES_H
#define OBLIVIOUSROUTING_ELECTRICAL_FLOW_PAR_BATCHES_H


#include <algorithm>
#include "../../datastructures/GraphADJ.h"
#include "../../utils/hash.h"
#include "../../solver/solver.h"
#include "../../solver/mwu_framework.h"
#include "../laplacian_solvers/AMGSolver.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "../../datastructures/IGraph.h"


class ElectricalFlowParallelBatches : public LinearObliviousSolverBase, public MWUFramework {

    int n, m;
    double K = 0; // this is used as an approximation errror for the Laplacian Solver( see. Paper for details)

    double epsilon = 0.5; // sketching parameter
    double epsilon_L = EPS; // Laplacian solving accuracy

    std::vector<std::unique_ptr<LaplacianSolver>> amg_pool;
    int p_threads = 1;

    double roh = 0.0;
    double alpha_local = 0.0;
    double inv_m = 0.0;
    std::vector<double> div_accum; // size n
    std::vector<int> tree_parent;

    std::vector<std::pair<int, int> > edges; // u<v only
    std::vector<double> edge_diffs;

    std::vector<double> edge_weights;             // w_e
    std::vector<double> edge_capacities;          // c_e
    std::vector<double> edge_distances;           // x_e
    std::vector<double> edge_probabilities;       // p_e


    double cap_X = 0.0;
    bool debug = false;
    int x_fixed = 0; // fixed node



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

    ElectricalFlowParallelBatches(IGraph& g, int root, bool debug = false):LinearObliviousSolverBase(g, root), n(g.getNumNodes()), m(g.getNumEdges()/2) {}

    void computeBasisFlows(LinearRoutingTable &table) override {
        init(debug);
        run(table);
        scaleFlowDown(table);
    }

    void run(LinearRoutingTable &table);
    void scaleFlowDown(LinearRoutingTable &table);

    void getApproxLoad(std::vector<double>& load);

    void init( bool debug = false, boost::property_tree::ptree _params = boost::property_tree::ptree());



};


#endif //OBLIVIOUSROUTING_ELECTRICAL_FLOW_PAR_BATCHES_H