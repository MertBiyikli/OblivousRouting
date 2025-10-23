//
// Created by Mert Biyikli on 06.08.25.
//

#ifndef OBLIVIOUSROUTING_AMGSOLVER_H
#define OBLIVIOUSROUTING_AMGSOLVER_H

#include <vector>
#include <tuple>
#include <iostream>
#include <unordered_set>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_builder.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
//#include <amgcl/adapter/crs.hpp>



#include "../../graph.h"
#include <unordered_map>
#include "../../utils/hash.h"
#include "LaplacianSolver.h"
#include <iostream>

class AMGSolver : public LaplacianSolver {
private:

    // AMG Solver define types
    typedef amgcl::make_solver<
            amgcl::amg<
                    amgcl::backend::builtin<double>,
                    amgcl::coarsening::smoothed_aggregation,
                    amgcl::relaxation::spai0
            >,
            amgcl::solver::cg<amgcl::backend::builtin<double>>
    > Solver;


    std::unique_ptr<Solver> solver;


public:

    void check_openmp_runtime();


    int n_reduced = 0;

 //   std::unordered_map<std::pair<int, int>, double> m_edge_weights; // Edge weights
 //   void init(std::unordered_map<std::pair<int, int>, double>& _edge_weights, int n, bool debug = false);
 //   void buildLaplacian();
    // void buildReducedLaplacian();
    std::vector<double> solve(const std::vector<double> &b) override;
    Eigen::VectorXd solve(const Eigen::VectorXd &b) override;
    // Eigen::VectorXd solve_reduced(const Eigen::VectorXd &b_full);
    void updateSolver() override;
    void buildLaplacian() override;
    void buildLaplacian_(const Graph& g) override;

    void updateAllEdges(const std::vector<double>& new_weights,
                        const std::vector<std::pair<int,int>>& edges) override {
        if (new_weights.size() != edges.size()) {
            throw std::runtime_error("updateAllEdges: size mismatch between weights and edges");
        }

        for (size_t e = 0; e < edges.size(); ++e) {
            int u = edges[e].first;
            int v = edges[e].second;
            double old_w = m_edge_weights[{u, v}];
            double new_w = new_weights[e];
            double delta = new_w - old_w;

            if (std::abs(delta) < 1e-18) continue;

            // Keep map symmetric
            m_edge_weights[{u, v}] = new_w;
            m_edge_weights[{v, u}] = new_w;

            // --- Update CSR Laplacian entries ---
            // Diagonal contributions
            m_values[m_indexMap[{u, u}]] += delta;
            m_values[m_indexMap[{v, v}]] += delta;

            // Off-diagonal contributions
            m_values[m_indexMap[{u, v}]] -= delta;
            m_values[m_indexMap[{v, u}]] -= delta;

        }
    }

    bool updateEdge(int u, int v, double new_weight);

    /*bool checkMatrix() const;*/

};

#endif //OBLIVIOUSROUTING_AMGSOLVER_H
