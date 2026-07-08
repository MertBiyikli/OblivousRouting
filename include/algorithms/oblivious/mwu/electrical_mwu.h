//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_ELECTRICAL_MWU_H
#define OBLIVIOUSROUTING_ELECTRICAL_MWU_H

#include "../oblivious_solver.h"
#include "mwu_framework.h"
#include "utils/time_tracking.h"
#include <Eigen/Sparse>
#include "oracle/electrical/laplacian_solver.h"
#include "core/errors.h"
/*
* This is the implementation of the electrical flow based MWU oblivious routing algorithm presented by Goranci et. al. in 2023.
* The idea is to repeatedly invoke an electrical flow computation (Laplacian solve)
* and then update the edge resistances based on the load on the edges. The process is repeated until convergence.
*/
class ElectricalMWU : public MWUFramework {
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

    bool use_sketching;

    std::vector<std::pair<int, int> > edges; // u<v only
    std::vector<double> edge_weights;             // w_e
    std::vector<double> edge_capacities;          // c_e
    std::vector<double> edge_distances;           // x_e
    std::vector<double> edge_probabilities;       // p_e
    std::vector<double> edge_diffs;

    // Preallocate as class members to avoid reallocs
    Eigen::SparseMatrix<double> X;          // n × ℓ precomputed RHS
    Eigen::MatrixXd SketchMatrix_t; // sketch matrix transposed

public:

    ElectricalMWU(IGraph& g, int root, bool use_sketching, bool debug = false)
    : MWUFramework(g, root), n(g.getNumNodes()), m(g.getNumUndirectedEdges()) {
        this->debug = debug;
        this->use_sketching = use_sketching;
    }

    // entry point
    virtual Result<void> computeBasisFlows(LinearRoutingTable& table) override {
        auto res = init(debug);
        if (!res) {
            return getError(res);
        }

        res = run(table);
        if (!res) {
            return getError(res);
        }

        res = scaleFlowDown(table);
        if (!res) {
            return getError(res);
        }

        return {};
    }

    virtual void printAdditionalStats() override {
        //nothing here..
    }

    virtual Result<void> updateDistances(const std::vector<double>& load) override;

    virtual Result<void> init(bool debug = false, boost::property_tree::ptree _params = boost::property_tree::ptree() );
    virtual Result<void> initAMGSolver(boost::property_tree::ptree _params);
    virtual Result<void> run(LinearRoutingTable &table);
    virtual Result<void> scaleFlowDown(LinearRoutingTable &table);

    virtual Result<void> getApproxLoad(std::vector<double>& load);
    Result<void> getExactLoad(std::vector<double>& load);



    // Helpers
    void extractEdges();
    void initEdgeDistances();
    Eigen::SparseMatrix<double> buildIncidence();
    Eigen::MatrixXd getSketchMatrix(double epsilon = 0.5);
    void addFlowToTable(const int& u, const Eigen::VectorXd& potential, LinearRoutingTable &table);
    void setEpsilon(double eps);

};


#endif //OBLIVIOUSROUTING_ELECTRICAL_MWU_H