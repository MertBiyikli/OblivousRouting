//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_LAPLACIAN_SOLVER_H
#define OBLIVIOUSROUTING_LAPLACIAN_SOLVER_H

#include <vector>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/coarsening/runtime.hpp>
#include <amgcl/relaxation/runtime.hpp>
#include "amgcl/coarsening/runtime.hpp"
#include "amgcl/solver/runtime.hpp"
#include <amgcl/make_solver.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <Eigen/Dense>

#include "graph_to_laplacian.h"
#include "../../../../utils/my_math.h"
#ifdef OR_ENABLE_OPENMP
#include <omp.h>
#endif

class LaplacianSolver {
protected:
    using Backend = amgcl::backend::builtin<double>;

    typedef amgcl::make_solver<
    amgcl::amg<
        Backend,
        amgcl::runtime::coarsening::wrapper,
        amgcl::runtime::relaxation::wrapper
        >,
    amgcl::runtime::solver::wrapper<Backend>
    > AMG;

    GraphToLaplacian weight_model;

    std::unique_ptr<AMG> hierarchy;
    bool use_dirichlet;
    int dirichlet_root;
    std::vector<double> m_values_dirichlet;

    bool debug;
    int n;
    // CSR storage
    std::vector<int> m_row_ptr;
    std::vector<int> m_col_ind;
    std::vector<double> m_values;

    // reusable buffers
    std::vector<double> result;
    std::vector<double> bvec_buffer, x_buffer;

    // this for configuration of the solver, e.g. coarsening and relaxation types
    boost::property_tree::ptree params;

public:
    LaplacianSolver() {
        n = 0;
        use_dirichlet = true;
        dirichlet_root = 0;
        debug = false;
    }
    virtual ~LaplacianSolver() = default;


    void init(IGraph& g, std::vector<double>& _adj_edge_weights, int n, const std::vector<std::pair<int, int>>& edges, bool debug = false);

    void updateSolver();
    void updateAllEdges(const std::vector<double> &new_weights, const std::vector<std::pair<int, int> > &edges);
    void buildLaplacian();

    std::vector<double> solve(const std::vector<double> &b, double eps = EPS);
    Eigen::VectorXd solve(const Eigen::VectorXd &b, double eps = EPS);


    void setSolverParams(const boost::property_tree::ptree& new_params);
    void print_params(const boost::property_tree::ptree& prm);

    void applyDirichletInPlace(std::vector<double>& vals);

    static void configureSingleThreadedAMG();
    static void printOpenMPDiagnostics(const std::string& context);
};

class OpenMPThreadGuard {
public:
    explicit OpenMPThreadGuard(int temporaryThreads) {
#ifdef OR_ENABLE_OPENMP
        previousThreads = omp_get_max_threads();
        previousDynamic = omp_get_dynamic();
        previousMaxActiveLevels = omp_get_max_active_levels();

        omp_set_dynamic(0);
        omp_set_max_active_levels(1);
        omp_set_num_threads(temporaryThreads);
#endif
    }

    ~OpenMPThreadGuard() {
#ifdef OR_ENABLE_OPENMP
        omp_set_dynamic(previousDynamic);
        omp_set_max_active_levels(previousMaxActiveLevels);
        omp_set_num_threads(previousThreads);
#endif
    }

private:
#ifdef OR_ENABLE_OPENMP
    int previousThreads = 1;
    int previousDynamic = 0;
    int previousMaxActiveLevels = 1;
#endif
};

#endif //OBLIVIOUSROUTING_LAPLACIAN_SOLVER_H