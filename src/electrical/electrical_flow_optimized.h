//
// Created by Mert Biyikli on 30.09.25.
//

#ifndef OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H
#define OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H

#include <algorithm>
#include "../graph.h"
#include "../utils/hash.h"
#include "../solver/solver.h"
#include "laplacian_solvers/AMGSolver.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "../tree_based/frt/raecke_frt_transform.h"
#include "laplacian_solvers/LaplacianSolverFactory.h"

struct WeightedEdge {
    std::pair<int, int> edge; // (u,v) with u<v
    double weight;
};

class ElectricalFlowOptimized : public ObliviousRoutingSolver{
    private:

    int n, m;
    // LaplacianSolver solver;
    Graph m_graph;
    std::unique_ptr<LaplacianSolver> amg;
    double roh = 0.0;
    double alpha_local = 0.0;
    double inv_m = 0.0;
    // std::unordered_map<std::pair<int, int>, double> x_edge2distance, p_edge2probability, w_edges2weights, c_edges2capacities;

    std::vector<double> edge_weights;             // w_e
    std::vector<double> edge_capacities;          // c_e
    std::vector<double> edge_distances;           // x_e
    std::vector<double> edge_probabilities;       // p_e

    // --- Per-adjacency (directed, aligned with Graph::neighbors(u)) ---
    std::vector<std::vector<int>> adj_edge_ids;   // maps (u, idx) -> edge_id (same id on both directions)


    //std::vector<std::vector<double>> x_edge_distance, p_edge_probability, w_edge_weight, c_edge_capacity;
    double cap_X = 0.0;

    // just for debugging purposes
    std::unordered_map<std::pair<int, int>, double> x_edge2distance, p_edge2probability, w_edges2weights, c_edges2capacities;


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

    std::vector<Eigen::SparseVector<double>>  getRoutingMatrix();
    Eigen::SparseMatrix<double> getSketchMatrix(int m, int n, double epsilon = 0.5);
    // double recoverNorm(const Eigen::MatrixXd& M, const Eigen::SparseVector<double>& vec);

    public:

    void solve(const Graph& g) override {
        this->init(g, debug);
        this->run();
    }
    void storeFlow() override;
    double getFlowForCommodity(int edge_id, int source, int target);

    // std::unordered_map<std::pair<int, int>, std::unordered_map<std::pair<int, int>, double>> f_e_st;
    std::unordered_map<std::pair<int, int>, double> f_e_u; // store the flow of the edge u→x for a fixed vertex x
    int x_fixed;

    // maps (u,v) → edge_id, stores both orientations
    inline long long make_key(int u, int v) const {
        if (u > v) std::swap(u,v); // always store with u<v
        return ( (static_cast<long long>(u) << 32) | static_cast<unsigned int>(v) );
    }

    std::unordered_map<long long, int> edge_id_map;

    // std::vector<WeightedEdge> edges_with_weights; // u<v only


    std::vector<std::pair<int, int> > edges; // u<v only
    // std::vector<double>              weights; // same order
    void run();
    void scaleFlowDown();

    void getApproxLoad(std::vector<double>& load);
    // Compute the median absolute difference between two node potentials
    // across all ℓ sampled solutions.
    //   diffs_u = potentials[u][*]  (length ℓ)
    //   diffs_v = potentials[v][*]  (length ℓ)
    double recoverNorm(const std::vector<double>& diffs_u,
                                            const std::vector<double>& diffs_v);

    void init(const Graph& g, bool debug = false, const std::string& solver_name = ("amg_cg"));



    // for debugging purposes
    bool CheckIfAdjacencyListInSyncWithUnordered_Map() const;

    std::vector<std::vector<int>> adj_f_e_u_id; // per-adjacency list version of f_e_u (stores edge ids)
    std::vector<std::vector<double>> adj_f_e_u; // per-adjacency list version of f_e_u
    void addFlow(int edge_id, int u, double flow) {
        // keep the edges and flows sorted by id
        if (adj_f_e_u_id[edge_id].empty() || adj_f_e_u_id[edge_id].back() < u) {
            adj_f_e_u_id[edge_id].push_back(u);
            adj_f_e_u[edge_id].push_back(flow);
        } else {
            auto it = std::ranges::lower_bound(adj_f_e_u_id[edge_id], u);
            int idx = std::distance(adj_f_e_u_id[edge_id].begin(), it);
            if (it != adj_f_e_u_id[edge_id].end() && *it == u) {
                adj_f_e_u[edge_id][idx] += flow; // accumulate
            } else {
                adj_f_e_u_id[edge_id].insert(it, u);
                adj_f_e_u[edge_id].insert(adj_f_e_u[edge_id].begin() + idx, flow);
            }
        }
    }

};

#endif //OBLIVIOUSROUTING_ELECTRICAL_FLOW_OPTIMIZED_H