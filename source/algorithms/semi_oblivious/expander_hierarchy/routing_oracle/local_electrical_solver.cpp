#include "algorithms/semi_oblivious/expander_hierarchy/routing_oracle/local_electrical_solver.h"


#include "core/errors.h"


#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>




#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

constexpr double kEmbeddingEpsilon = 1e-8;

// TODO: replace this with the Laplacian Solver
Result<std::vector<double>> solveGroundedLaplacian(const int node_count,const std::vector<std::pair<int, int>>& edges,const std::vector<double>& conductances,const std::vector<double>& imbalance,const int grounded_vertex) {
    if (node_count <= 0) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Grounded Laplacian requires at least one vertex.");
    }

    if (grounded_vertex < 0 ||grounded_vertex >= node_count) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Grounded vertex is outside the local graph.");
    }

    if (edges.size() != conductances.size()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Edge and conductance vectors have different sizes.");
    }

    if (imbalance.size() != static_cast<std::size_t>(node_count)) {
        return makeErrorMessage(ErrorCode::InvalidDemand,"Imbalance vector has an invalid size.");
    }

    const double imbalance_sum = std::accumulate(imbalance.begin(),imbalance.end(),0.0);

    if (!std::isfinite(imbalance_sum) || std::abs(imbalance_sum) > 1e-9) {
        return makeErrorMessage(ErrorCode::InvalidDemand,"Electrical demand must be balanced. Sum=" +std::to_string(imbalance_sum));
    }

    for (const double value : imbalance) {
        if (!std::isfinite(value)) {
            return makeErrorMessage(ErrorCode::InvalidDemand,"Imbalance vector contains a non-finite value.");
        }
    }

    /*
     * A one-vertex graph has only the grounded potential.
     * Its balanced demand must be zero.
     */
    if (node_count == 1) {
        return std::vector<double>{0.0};
    }

    /*
     * Remove the grounded vertex from the linear system.
     *
     * reduced_index[v] gives the row/column of v in the
     * reduced (n - 1) x (n - 1) Laplacian.
     */
    std::vector<int> reduced_index(node_count, -1);

    int next_reduced_index = 0;

    for (int v = 0;v < node_count;++v) {
        if (v == grounded_vertex) {
            continue;
        }
        reduced_index[v] =next_reduced_index++;
    }

    const int reduced_node_count =node_count - 1;

    std::vector<double> diagonal(reduced_node_count,0.0);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(2 * edges.size() +static_cast<std::size_t>(reduced_node_count));

    /*
     * Construct the grounded Laplacian.
     *
     * For an edge (u,v) with conductance c:
     *
     *     L[u,u] += c
     *     L[v,v] += c
     *     L[u,v] -= c
     *     L[v,u] -= c
     *
     * Entries involving the grounded vertex are omitted from
     * the reduced matrix, but their conductance still contributes
     * to the diagonal of the non-grounded endpoint.
     */
    for (std::size_t e = 0;e < edges.size();++e) {
        const auto [u, v] = edges[e];
        const double conductance = conductances[e];

        if (u < 0 || v < 0 ||u >= node_count ||v >= node_count) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Laplacian edge contains an invalid endpoint.");
        }

        if (u == v) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Laplacian graph contains a self-loop.");
        }

        if (!std::isfinite(conductance) ||conductance <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Laplacian edge contains a non-positive conductance.");
        }

        if (u != grounded_vertex) {
            diagonal[reduced_index[u]] +=conductance;
        }

        if (v != grounded_vertex) {
            diagonal[reduced_index[v]] +=conductance;
        }

        if (u != grounded_vertex &&
            v != grounded_vertex) {
            triplets.emplace_back(
                reduced_index[u],
                reduced_index[v],
                -conductance
            );

            triplets.emplace_back(
                reduced_index[v],
                reduced_index[u],
                -conductance
            );
        }
    }

    for (int row = 0;row < reduced_node_count;++row) {
        if (!std::isfinite(diagonal[row]) ||diagonal[row] <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Grounded Laplacian contains a non-positive diagonal. The graph may be disconnected."
            );
        }

        triplets.emplace_back(row,row,diagonal[row]);
    }

    Eigen::SparseMatrix<double> grounded_laplacian(reduced_node_count,reduced_node_count);

    grounded_laplacian.setFromTriplets(triplets.begin(),triplets.end());

    grounded_laplacian.makeCompressed();

    Eigen::VectorXd reduced_rhs(reduced_node_count);

    for (int v = 0;v < node_count;++v) {
        if (v == grounded_vertex) {
            continue;
        }

        reduced_rhs[reduced_index[v]] =imbalance[v];
    }



    /*
     * A grounded Laplacian is symmetric positive definite when
     * the graph is connected, so SimplicialLDLT is appropriate.
     */
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    solver.compute(grounded_laplacian);

    if (solver.info() != Eigen::Success) {
        return makeErrorMessage(ErrorCode::SolverFailed,"Could not factorize the grounded Laplacian. The local graph may be disconnected.");
    }

    const Eigen::VectorXd reduced_potential = solver.solve(reduced_rhs);

    if (solver.info() != Eigen::Success) {
        return makeErrorMessage(ErrorCode::SolverFailed,"Grounded Laplacian solve failed.");
    }

    if (!reduced_potential.allFinite()) {
        return makeErrorMessage(ErrorCode::SolverFailed,"Grounded Laplacian returned non-finite potentials.");
    }

    /*
     * Restore the full potential vector.
     *
     * Potentials are unique only up to an additive constant.
     * We select the representative with:
     *
     *     potential[grounded_vertex] = 0.
     */
    std::vector<double> potential(node_count,0.0);

    for (int v = 0;v < node_count;++v) {
        if (v == grounded_vertex) {
            potential[v] = 0.0;
            continue;
        }

        potential[v] = reduced_potential[reduced_index[v]];
    }

    return potential;
}


Result<std::vector<double>> solveGroundedLaplacianViaAMGCL(const int node_count,const std::vector<std::pair<int, int>>& edges,const std::vector<double>& conductances,const std::vector<double>& imbalance,const int grounded_vertex) {
    if (node_count <= 0) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Grounded Laplacian requires at least one vertex.");
    }

    if (grounded_vertex < 0 ||grounded_vertex >= node_count) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Grounded vertex is outside the local graph.");
    }

    if (edges.size() != conductances.size()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Edge and conductance vectors have different sizes.");
    }

    if (imbalance.size() != static_cast<std::size_t>(node_count)) {
        return makeErrorMessage(ErrorCode::InvalidDemand,"Imbalance vector has an invalid size.");
    }

    const double imbalance_sum = std::accumulate(imbalance.begin(),imbalance.end(),0.0);

    if (!std::isfinite(imbalance_sum) || std::abs(imbalance_sum) > 1e-9) {
        return makeErrorMessage(ErrorCode::InvalidDemand,"Electrical demand must be balanced. Sum=" +std::to_string(imbalance_sum));
    }

    for (const double value : imbalance) {
        if (!std::isfinite(value)) {
            return makeErrorMessage(ErrorCode::InvalidDemand,"Imbalance vector contains a non-finite value.");
        }
    }

    /*
     * A one-vertex graph has only the grounded potential.
     * Its balanced demand must be zero.
     */
    if (node_count == 1) {
        return std::vector<double>{0.0};
    }

    /*
     * Remove the grounded vertex from the linear system.
     *
     * reduced_index[v] gives the row/column of v in the
     * reduced (n - 1) x (n - 1) Laplacian.
     */
    std::vector<int> reduced_index(node_count, -1);

    int next_reduced_index = 0;

    for (int v = 0;v < node_count;++v) {
        if (v == grounded_vertex) {
            continue;
        }
        reduced_index[v] =next_reduced_index++;
    }

    const int reduced_node_count =node_count - 1;

    std::vector<double> diagonal(reduced_node_count,0.0);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(2 * edges.size() +static_cast<std::size_t>(reduced_node_count));

    LaplacianSolver amgcl;
    auto initialize = amgcl.initGrounded(node_count, edges, conductances, grounded_vertex);
    if (!initialize) {
        return getError(initialize);
    }
    auto pot = amgcl.solve(imbalance);

    if (!pot) {
        return getError(pot);
    }

    if (pot->size() != static_cast<std::size_t>(node_count)) {
        return makeErrorMessage(
            ErrorCode::SolverFailed,
            "AMGCL returned an invalid potential-vector size."
        );
    }
    /*
     * Enforce the selected representative explicitly.
     */
    (*pot)[grounded_vertex] = 0.0;

    return pot;
}

Result<LocalElectricalFlowResult> LocalElectricalSolver::routeDemand(const IGraph& graph,const HierarchyCluster& cluster,const std::vector<double>& local_imbalance) const {
    const int local_node_count = static_cast<int>(cluster.original_vertices.size());

    if (local_node_count == 0) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Cannot solve electrical flow on an empty cluster.");
    }

    if (local_imbalance.size() != cluster.original_vertices.size()) {
        return makeErrorMessage(ErrorCode::InvalidDemand,"Cluster imbalance has invalid size.");
    }

    const double imbalance_sum = std::accumulate(local_imbalance.begin(),local_imbalance.end(),0.0);

    if (!std::isfinite(imbalance_sum) || std::abs(imbalance_sum) > kEmbeddingEpsilon) {
        return makeErrorMessage(ErrorCode::InvalidDemand,"Cluster electrical demand is not balanced. Sum=" +std::to_string(imbalance_sum));
    }

    bool trivial = true;

    for (const double value : local_imbalance) {
        if (!std::isfinite(value)) {
            return makeErrorMessage(ErrorCode::InvalidDemand,"Cluster demand contains a non-finite value.");
        }

        if (std::abs(value) > kEmbeddingEpsilon) {
            trivial = false;
        }
    }

    if (trivial) {
        return LocalElectricalFlowResult{
            .edge_flow =std::vector<double>(cluster.induced_edges.size(),0.0),
            .max_congestion = 0.0,
            .max_conservation_error = 0.0
        };
    }

    std::unordered_map<int, int> global_to_local;
    global_to_local.reserve(cluster.original_vertices.size());

    for (int local = 0;local < local_node_count;++local) {
        global_to_local.emplace(cluster.original_vertices[local],local);
    }

    std::vector<std::pair<int, int>> local_edges;
    std::vector<double> conductances;
    std::vector<int> local_to_cluster_edge;

    local_edges.reserve(cluster.induced_edges.size());
    conductances.reserve(cluster.induced_edges.size());
    local_to_cluster_edge.reserve(cluster.induced_edges.size());

    for (int cluster_edge_index = 0; cluster_edge_index <static_cast<int>(cluster.induced_edges.size());++cluster_edge_index) {
        const int global_edge_id = cluster.induced_edges[cluster_edge_index];

        if (global_edge_id < 0 || global_edge_id >= graph.getNumDirectedEdges()) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Cluster contains an invalid global edge ID.");
        }

        const auto [global_u, global_v] = graph.getEdgeEndpoints(global_edge_id);

        const auto u_it = global_to_local.find(global_u);
        const auto v_it = global_to_local.find(global_v);

        if (u_it == global_to_local.end() ||v_it == global_to_local.end()) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Cluster induced edge leaves the cluster.");
        }

        const double conductance = graph.getEdgeCapacity(global_edge_id);

        if (!std::isfinite(conductance) ||conductance <= 0.0) {
            return makeErrorMessage(ErrorCode::InvalidGraph,"Cluster edge has invalid conductance.");
        }

        local_edges.emplace_back(
            u_it->second,
            v_it->second
        );

        conductances.push_back(conductance);
        local_to_cluster_edge.push_back(cluster_edge_index);
    }

    if (local_edges.empty()) {
        return makeErrorMessage(ErrorCode::InvalidGraph,"Nontrivial cluster demand has no internal edges.");
    }

    /*
     * Any local vertex may be grounded. Grounding removes the
     * one-dimensional nullspace of the graph Laplacian.
     */
    const int grounded_vertex = 0;

    auto potential_result = solveGroundedLaplacian(
        local_node_count,
        local_edges,
        conductances,
        local_imbalance,
        grounded_vertex
    );

    if (!potential_result) {
        return getError(potential_result);
    }



    const auto& potential = *potential_result;

    if (potential.size() !=
        cluster.original_vertices.size()) {
        return makeErrorMessage(ErrorCode::SolverFailed,"Electrical solver returned invalid potential size.");
    }

    LocalElectricalFlowResult result;
    result.edge_flow.assign(cluster.induced_edges.size(),0.0);

    std::vector<double> actual_divergence(local_node_count,0.0);

    for (int local_edge_id = 0;local_edge_id <static_cast<int>(local_edges.size()); ++local_edge_id) {
        const auto [u, v] =local_edges[local_edge_id];
        const double conductance =conductances[local_edge_id];

        const double flow = conductance *(potential[u] - potential[v]);

        if (!std::isfinite(potential[u]) || !std::isfinite(potential[v])) {
            return makeErrorMessage(
                ErrorCode::SolverFailed,
                "AMGCL returned non-finite potentials in cluster " +
                    std::to_string(cluster.id) +
                    ". Local edge=" +
                    std::to_string(local_edge_id) +
                    ", endpoints=(" +
                    std::to_string(u) +
                    "," +
                    std::to_string(v) +
                    "), potential_u=" +
                    std::to_string(potential[u]) +
                    ", potential_v=" +
                    std::to_string(potential[v])
            );
    }

        const int cluster_edge_index = local_to_cluster_edge[local_edge_id];

        result.edge_flow[cluster_edge_index] = flow;

        actual_divergence[u] += flow;
        actual_divergence[v] -= flow;

        result.max_congestion = std::max(
            result.max_congestion,
            std::abs(flow) / conductance
        );
    }

    for (int local_vertex = 0;local_vertex < local_node_count;++local_vertex) {
        result.max_conservation_error = std::max(result.max_conservation_error,std::abs(actual_divergence[local_vertex] -local_imbalance[local_vertex])
        );
    }

    if (result.max_conservation_error > kEmbeddingEpsilon) {
        return makeErrorMessage(ErrorCode::SolverFailed,"Local electrical conservation failed in cluster " +std::to_string(cluster.id) +" with error " +std::to_string(result.max_conservation_error));
    }

    return result;
}