//
// Created by Mert Biyikli on 09.07.25.
//



#include "electrical_flow_naive.h"
#include <random>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

// Checks for issues in a CSR Laplacian matrix representation
void diagnose_csr_laplacian(
        const std::vector<int>& row_ptr,
        const std::vector<int>& col_ind,
        const std::vector<double>& values,
        double tolerance = 1e-6
) {
    int n = static_cast<int>(row_ptr.size()) - 1;

    std::vector<double> row_sums(n, 0.0);
    std::vector<bool> nonzero_diag(n, false);
    std::vector<bool> nonzero_row(n, false);
    std::vector<bool> nonzero_col(n, false);

    for (int i = 0; i < n; ++i) {
        for (int idx = row_ptr[i]; idx < row_ptr[i + 1]; ++idx) {
            int j = col_ind[idx];
            double v = values[idx];

            row_sums[i] += v;

            if (i == j && std::abs(v) > tolerance)
                nonzero_diag[i] = true;

            if (std::abs(v) > tolerance) {
                nonzero_row[i] = true;
                nonzero_col[j] = true;
            }
        }
    }

    bool issues = false;
    for (int i = 0; i < n; ++i) {
        if (!nonzero_diag[i]) {
            std::cerr << "⚠️  Zero or missing diagonal at row " << i << std::endl;
            issues = true;
        }
        if (!nonzero_row[i]) {
            std::cerr << "⚠️  Entire zero row at " << i << std::endl;
            issues = true;
        }
        if (!nonzero_col[i]) {
            std::cerr << "⚠️  Entire zero column at " << i << std::endl;
            issues = true;
        }
        if (std::abs(row_sums[i]) > tolerance) {
            std::cerr << "⚠️  Row " << i << " has non-zero sum: " << row_sums[i] << std::endl;
            issues = true;
        }
    }

    if (!issues) {
        std::cout << "✅ Matrix passed all CSR diagnostics. Looks valid." << std::endl;
    }
}

double NegativeExponent(double base, int exp) {
    if (base == 0.0) {
        throw std::invalid_argument("Base cannot be zero for negative exponent.");
    }
    if (exp < 0) {
        throw std::invalid_argument("Exponent must be non-negative.");
    }
    return 1.0 / std::pow(base, exp);
}




void ElectricalFlowNaive::init(const Graph &g, bool debug) {
    m_graph = g;

    n = m_graph.getNumNodes();
    m = m_graph.getNumEdges();

    // set algorithm parameters
    roh = std::sqrt(2.0*static_cast<double>(m)); // Initialize roh based on the number of nodes
    alpha_local = std::log2(n)*std::log2(n); // Initialize alpha_local based on the number of nodes
    this->cap_X = m;
    this->number_of_iterations = 8.0*roh*std::log(m)/alpha_local;
    this->inv_m = 1.0 / static_cast<double>(m);


    // fix a node x
    this->x_fixed = rand()%n; // Randomly select a fixed node x from the graph

    // initialize distance, weights, probabilities
    initEdgeDistances();
    extract_edge_list();

    // initialize the Laplacian Solver
    amg.init(w_edges2weights, n, debug);
    this->debug = debug;

    // compute diagonal matrix consisting of only the capacities
    U = Eigen::SparseMatrix<double>(m, m);
    for(int i = 0; i < m; ++i) {
        U.insert(i, i) = c_edges2capacities[edges[i]]; // Insert the capacities into the diagonal matrix
    }
}

void ElectricalFlowNaive::extract_edge_list() {
    edges.clear();
    weights.clear();
    for(auto &kv : x_edge2distance) {
        auto [u,v] = kv.first;
        if(u < v) {
            edges.emplace_back(u,v);
            weights.push_back(kv.second);
        }
    }

    // sort the edges based on the first node, then second node
    std::sort(edges.begin(), edges.end(),
              [](const std::pair<int, int> &a, const std::pair<int, int> &b) {
                  if (a.first != b.first) return a.first < b.first;
                  return a.second < b.second;
              });
}


void ElectricalFlowNaive::initEdgeDistances() {
    // Initialize the edge distances and probabilities
    for(int u = 0; u < n; ++u) {
        for(int v : m_graph.neighbors(u)) {
            double distance = 1.0;
            x_edge2distance[{u, v}] = distance; // Initialize distances to 1.0
            c_edges2capacities[{u, v}] = m_graph.getEdgeCapacity( u, v);

            // Update the probability based on the new distance
            double p = distance/cap_X;
            p_edge2probability[{u, v}] = p;

            double w = std::pow(c_edges2capacities[{u, v}], 2) * 1.0/(p + inv_m); // Calculate the weight based on the probability and inverse of the number of edges
            w_edges2weights[{u, v}] = w;
        }
    }

}


void ElectricalFlowNaive::run() {

    Eigen::VectorXd demand_u_x = Eigen::VectorXd::Zero(m_graph.getNumNodes()); // Initialize the demand vector for the fixed node x

    for(int t = 0; t < std::min(number_of_iterations, std::numeric_limits<int>::max()); ++t) {

        // TODO: avoid using matrix multiplication for the routing matrix whatsever
        // since you can avoid this by storing the flow values directly as an map
        // 3) Build M = W B (Bᵀ W B)^{-1} via m computations
        Eigen::SparseMatrix<double> M = getRoutingMatrix();

        // TODO: shift this to an extra function that is called within getRoutingMatrix()
        // update the flow values based on new routing matrix
        for(int u = 0; u < n; ++u) {
            if(u == x_fixed) continue; // Skip the fixed node x

            demand_u_x[u] = 1.0;
            demand_u_x[x_fixed] = -1.0;
            Eigen::VectorXd flow = M*demand_u_x; // Get the flow for the edge


            for(int flow_edge_id = 0; flow_edge_id < flow.size(); ++flow_edge_id) {
                double flow_value = flow[flow_edge_id];
                if(std::abs(flow_value) > 1e-16) { // Only update if the flow is significant
                    f_e_u[{flow_edge_id, u}] += flow_value; // Update the flow for the edge u→x
                }
            }
            demand_u_x[u]=0.0;
        }

        // print the flow values for debugging
        if(debug) {
            for(auto& [edge, flow] : f_e_u) {
                std::cout << "Flow for edge (" << edges[edge.first].first << ", " << edges[edge.first].second << ") commodity ( " << edge.second << " -> " << x_fixed << " ):" << flow << "\n";
            }
        }

        Ms.push_back(M); // Store the M matrix for later use

        if (debug){
            std::cout << "M has " << M.nonZeros() << " nonzeros\n";

            // print the M matrix for debugging
            std::cout << "M matrix (first 5 rows and columns):\n";
            int n = M.rows(), m = M.cols();
            for (int i = 0; i < std::min(5, n); ++i) {
                for (int j = 0; j < std::min(5, m); ++j) {
                    std::cout << M.coeff(i, j) << " ";
                }
                std::cout << "\n";
            }
        }

        auto aload = getApproxLoad();
        // print approxiload
        if(debug) {
            std::cout << "Approximate load:\n";
            for (const auto &kv: aload) {
                auto [edge, load] = kv;
                std::cout << "Edge (" << edge.first << ", " << edge.second << "): Load = " << load << "\n";
            }
        }

        cap_X = 0.0;
        // update the edge distances based on the approximate load
        for(const auto& [edge, load] : aload) {

            if(load <= 0.0 || std::isnan(load)) continue; // Skip invalid loads

            if(debug) {
                std::cout << "Updating edge (" << edge.first << ", " << edge.second << ") with load = " << load << "\n";
                std::cout << "Current distance: " << x_edge2distance[edge] << "\n";
            }
            x_edge2distance[edge] *= (1+(1/(2*roh))*load); // Update the distance with the approximate load
            cap_X += x_edge2distance[edge]; // Update the total capacity
            if(debug) {
                std::cout << "New distance: " << x_edge2distance[edge] << "\n";
            }
        }
        if(debug) {
            std::cout << "Total capacity X: " << cap_X << "\n";
        }

        updateEdgeDistances();

        // break;
    }


    // take the average of the Ms matrices
    M_avg = Ms[0];
    for(int i = 1; i < Ms.size(); ++i) {
        M_avg += Ms[i];
    }
    M_avg /= static_cast<double>(Ms.size());


    //take the average of the flows
    for(auto& [edge, flow] : f_e_u) {
        flow /= static_cast<double>(number_of_iterations); // Average the flow values
    }





    if(debug) {
        // print the average M matrix
        std::cout << "Average M matrix (first 5 rows and columns):\n";
        int n = M_avg.rows(), m = M_avg.cols();
        for (int i = 0; i < std::min(5, n); ++i) {
            for (int j = 0; j < std::min(5, m); ++j) {

                // only print out non-zero values, or values that are not too close to zero
                // to avoid printing too many zeros
                if(std::abs(M_avg.coeff(i, j))-1e-16<0) {
                    std::cout << "0 ";
                }else {
                    std::cout << M_avg.coeff(i, j) << " ";
                }
            }
            std::cout << "\n";
        }
        // print the average flow values
        std::cout << "Average flow values (first 5 edges):\n";
        for (auto &[e_u, flow]: f_e_u) {
            if (flow != 0.0) { // Only print significant flows
                std::cout << "Edge (" << edges[e_u.first].first << ", " << edges[e_u.first].second << ") commodity ( "
                          << e_u.second << " -> " << x_fixed << " ): Flow = " << flow << "\n";
            }
        }
    }



    // TODO: move this into an extra function that is called after the run() method
    // scale the flow values to meet one unit of flow per commodity
    // for this sum up the outgoing for each source s of each commodity pair (s, t)
    // and multiply \forall e \in E f_st(e) by 1/sum_{s} f_st(e)
    for (int s = 0; s<n; ++s) {
        if(s == x_fixed) continue; // Skip the fixed node x
        int source = s; // Get the source node of the commodity

        double outgoing_flow = 0.0; // Get the outgoing flow for the edge
        for(auto& neig : m_graph.neighbors(source)) {
            // if(source > neig) continue; // Ensure we only consider one direction of the edge

            auto edge = std::find(edges.begin(), edges.end(),
                                  (source < neig ?
                                    std::make_pair(source, neig) :
                                    std::make_pair(neig, source)));

            // Find the edge in the edges vector
            if(edge != edges.end()) {
                int edge_id = std::distance(edges.begin(), edge);
                outgoing_flow += std::abs(f_e_u[{edge_id, source}]); // Sum the flow for the edge
            }
        }

        // scale for each edge the flow by 1/outgoing_flow
        if (std::abs(outgoing_flow) > 1e-16) { //
            for(auto& [edge, flow] : f_e_u) {
                if(edge.second == source) { // Only scale the flow for the edges outgoing from the source
                    flow /= outgoing_flow; // Scale the flow by the outgoing flow
                }
            }
        } else {
            std::cerr << "Warning: No outgoing flow for source " << source << ". Skipping scaling.\n";
        }
    }

    // print the scaled flow values
    if(debug) {
        std::cout << "Scaled flow values (first 5 edges):\n";
        for (auto &[e_u, flow]: f_e_u) {
            if (flow != 0.0) { // Only print significant flows
                std::cout << "Edge (" << edges[e_u.first].first << ", " << edges[e_u.first].second << ") commodity ( "
                          << e_u.second << " -> " << x_fixed << " ): Flow = " << flow << "\n";
            }
        }
    }
}


std::unordered_map<std::pair<int, int>, double> ElectricalFlowNaive::getApproxLoad() {

    std::unordered_map<std::pair<int, int>, double> approx_load;
    // build the incidence matrix B and the weight vector w
    if(!B.nonZeros()) {
        buildIncidence();
    }


    if (debug)
        std::cout << "Running sketch matrix generation...\n";
    auto C = getSketchMatrix(m, n, 0.5);


    auto X = B.transpose() * U * C.transpose();


    // TODO: could we store V as a Sparse Matrix instead??
    Eigen::MatrixXd V(n, X.cols() );

    // for each column of X, solve the linear system with the Laplacian of the graph itself
    for(int i = 0; i<X.cols(); ++i) {
        Eigen::VectorXd x = X.col(i);

        // TODO remove here the solver
        Eigen::VectorXd u = amg.solve(x);
        V.col(i) = u;
        if(debug) {
            std::cout << "rhs: \n" << x.transpose() << "\n";
            std::cout << "Column " << i << " of U:\n" << u.transpose() << "\n";
        }
    }

    // for each edge (u,v), compute the approximate load
    for(int i = 0; i < edges.size(); ++i) {
        auto b = B.row(i);

        double norm = recoverNorm(V.transpose(), b);

        approx_load[edges[i]] = w_edges2weights[edges[i]] * norm*(1/c_edges2capacities[edges[i]]); // Approximate load is the weight times the norm divided by the capacity

        if(debug) {
            std::cout << "Edge (" << edges[i].first << ", " << edges[i].second << "): "
                      << "Approximate Load = " << approx_load[edges[i]] << ", "
                      << "Weight = " << w_edges2weights[edges[i]] << ", "
                      << "Norm = " << norm << ", "
                      << "Capacity = " << c_edges2capacities[edges[i]] << "\n";
        }
    }


    return approx_load;
    // get the approxi load from the recover norm
}


Eigen::SparseMatrix<double>  ElectricalFlowNaive::getSketchMatrix(int _m, int _n, double _epsilon) {
    // ℓ = O(ε⁻² log(1/δ))
    double c = 1.1;
    double delta = NegativeExponent(_n, 10); // Set delta to a small value or a default value
    double epsilon = _epsilon; // Set epsilon to the provided value or a default value
    if(debug) {
        std::cout << "Generating sketch matrix with parameters:\n"
                  << "c = " << c << ", delta = " << delta << ", epsilon = " << epsilon
                  << ", m = " << _m << ", n = " << _n << "\n";
    }

    int l = ((c/(epsilon*epsilon) * std::log(1.0/delta) ));

    if(debug) {
        std::cout << "Generating sketch matrix with l = " << l << ", m = " << _m << ", n = " << _n << ", epsilon = " << epsilon << "\n";
    }

    Eigen::SparseMatrix<double> C(l, _m);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::cauchy_distribution<double> dist(0.0, 1.0);

    for(int i = 0; i < l; ++i) {
        for (int j = 0; j < _m; ++j) {
            double rand_value = dist(gen);
            // assign sparse matrix C with random values from the Cauchy distribution
            C.insert(i, j) = rand_value;
        }
    }

    // print sketch matrix for debugging
    if(debug) {
        std::cout << "Sketch matrix C (first 5 rows and columns):\n";
        int rows = C.rows(), cols = C.cols();
        for (int i = 0; i < std::min(5, rows); ++i) {
            for (int j = 0; j < std::min(5, cols); ++j) {
                std::cout << C.coeff(i, j) << " ";
            }
            std::cout << "\n";
        }
    }

    return C;
}


// compute the median of the Matrix M vector vec product
double ElectricalFlowNaive::recoverNorm(const Eigen::MatrixXd& M, const Eigen::VectorXd& vec) {
    auto s = M*vec;

    std::vector<double> abs_vals(s.size());
    for (int i = 0; i < s.size(); ++i) {
        abs_vals[i] = std::abs(s[i]);
    }
    size_t mid = abs_vals.size() / 2;
    std::nth_element(abs_vals.begin(), abs_vals.begin() + mid, abs_vals.end());
    return abs_vals[mid];
}


void ElectricalFlowNaive::updateEdgeDistances() {
    // Update the edge distances and probabilities
    for(auto& edge : edges) {
        int v = edge.first;
        int u = edge.second;

            double distance = x_edge2distance[edge];
            if(distance <= 0.0 || std::isnan(distance)) continue;

            // Update the probability based on the new distance
            double p = distance/cap_X;
            p_edge2probability[edge] = p;

            double w = std::pow(c_edges2capacities[edge], 2) * 1.0/(p + inv_m); // Calculate the weight based on the probability and inverse of the number of edges
            w_edges2weights[edge] = w;

            w_edges2weights[{u, v}] = w; // Update the weight for the edge in the same order as edges vector

            //amg.updateEdge(v, u, w);
            //solver.NonSyncUpdateEdgeWeights({v, u}, w); // Update the weights in the Laplacian solver

    }
    //solver.SyncEdgeWeights(); // Synchronize the weights in the Laplacian solver

    // TODO: change this by only changing the numeric values of the matrix solver
    amg.init(w_edges2weights, n);
    // amg.updateSolver();

    // Refresh our weights[] array in the same edge‐order
    for(int i = 0; i < edges.size(); ++i) {
        auto e = edges[i];
        weights[i] = w_edges2weights[e];
    }
    refreshWeightMatrix();


    if(debug) {
        std::cout << "Updated edge distances and weights.\n";
        std::cout << "Total capacity cap_X: " << cap_X << "\n";
        std::cout << "Number of edges: " << edges.size() << "\n";
        for(int i = 0; i < edges.size(); ++i) {
            auto e = edges[i];
            std::cout << "Edge (" << e.first << ", " << e.second << "): "
                      << "Distance = " << x_edge2distance[e] << ", "
                      << "Weight = " << w_edges2weights[e] << ", "
                      << "Probability = " << p_edge2probability[e] << ", "
                      << "Capacity = " << c_edges2capacities[e] << "\n";
        }
    }
}

void ElectricalFlowNaive::refreshWeightMatrix() {
    // Rebuild the weight diagonal matrix W
    if(!W.nonZeros()) {
        buildWeightDiag();
    } else {
        // If W already exists, just update it with the new weights
        for (int e = 0; e < m; ++e) {
            auto edge = edges[e];
            double weight = w_edges2weights[edge];
            W.coeffRef(e, e) = weight; // Update the diagonal entry
        }
    }

    if(debug) {
        std::cout << "Weight diagonal matrix W (sparse):\n";
        std::cout << "W has " << W.nonZeros() << " nonzeros\n";
    }
}




Eigen::SparseMatrix<double> ElectricalFlowNaive::getRoutingMatrix() {
    buildIncidence();
    buildWeightDiag();

    auto bw = W*B;
    auto bw_T = bw.transpose();

    Eigen::MatrixXd M(n, m);
    Eigen::MatrixXd L_inv(n, n);


    if(debug)
        diagnose_csr_laplacian(amg.m_row_ptr, amg.m_col_ind, amg.m_values);

    for(int i = 0; i<m; i++) {

        if(debug) {
            // print out the amg solver
            for(auto [e, w] : amg.m_edge_weights) {
                std::cout << "Edge: " << e.first << " " << e.second << " : " << w << std::endl;
            }

            std::cout << "row_ptr "<< std::endl;
            for(auto U : amg.m_row_ptr) {
                std::cout << "" << U << " ";
            }
            std::cout << std::endl;

            std::cout << "m_col_ind "<< std::endl;
            for(auto U : amg.m_col_ind) {
                std::cout << "" << U << " ";
            }
            std::cout << std::endl;

            std::cout << "m_values "<< std::endl;
            for(auto U : amg.m_values) {
                std::cout << "" << U << " ";
            }
            std::cout << std::endl;

        }

        Eigen::VectorXd rhs = bw_T.col(i);
        rhs.array() -= rhs.mean();

        if(debug) {
            std::cout << "Solving for column " << i << " with rhs: " << rhs.transpose() << std::endl;
        }

        auto result = amg.solve(rhs);
        M.col(i) = result;

        if(debug) {

            std::cout << "result for iteration: " << i << "\n";
            for (auto r: result) std::cout << r << " ";
            std::cout << std::endl;

            // print the Laplacian for debugging
            Eigen::SparseMatrix<double> Lsp(n, n);
            std::vector<Eigen::Triplet<double>> triplets;
            for (int row = 0; row < n; ++row) {
                for (int idx = amg.m_row_ptr[row]; idx < amg.m_row_ptr[row + 1]; ++idx) {
                    int col = amg.m_col_ind[idx];
                    double val = amg.m_values[idx];
                    triplets.emplace_back(row, col, val);
                }
            }
            Lsp.setFromTriplets(triplets.begin(), triplets.end());
            std::cout << "Laplacian matrix L (sparse, first 5 rows and columns):\n";
            int rows = Lsp.rows(), cols = Lsp.cols();
            for (int r = 0; r < std::min(5, rows); ++r) {
                for (int c = 0; c < std::min(5, cols); ++c) {
                    std::cout << Lsp.coeff(r, c) << " ";
                }
                std::cout << "\n";
            }
        }
    }

    return M.transpose().sparseView();

}


// Build the oriented incidence matrix B (E×V).  Each edge e=(u,v):
//   row e has +1 in column u, –1 in column v.
void ElectricalFlowNaive::buildIncidence() {
    if(B.nonZeros() == 0) {
        B = Eigen::SparseMatrix<double>(m, n);
        std::vector<Eigen::Triplet<double>> T;
        T.reserve(2*m);
        for(int e = 0; e < m; ++e){
            auto [u,v] = edges[e];
            T.emplace_back(e, u, -1.0);
            T.emplace_back(e, v, +1.0);
        }
        B.setFromTriplets(T.begin(), T.end());
    }

    if(debug) {
        // Print out the edges
        int id = 0;
        for(auto& edge : edges) {
            std::cout << "Edge: " << (id++) << " (" << edge.first << ", " << edge.second << ")\n";
        }

        // print the incidence matrix B for debugging
        std::cout << "Incidence matrix B (first 5 rows and columns):\n";
        int rows = B.rows(), cols = B.cols();
        for (int i = 0; i < std::min(5, rows); ++i) {
            for (int j = 0; j < std::min(5, cols); ++j) {
                std::cout << B.coeff(i, j) << " ";
            }
            std::cout << "\n";
        }
    }
}




void ElectricalFlowNaive::buildWeightDiag() {
    if(!W.nonZeros()) {
        // 1. Build W as a sparse diagonal matrix
        W = Eigen::SparseMatrix<double>(m, m);
        std::vector<Eigen::Triplet<double>> w_triplets(m);
        for (int e = 0; e < m; ++e) {
            auto edge = edges[e];
            double weight = w_edges2weights[edge];
            w_triplets.push_back({e, e, weight}); // Add a triplet for the diagonal entry
        }
        W.setFromTriplets(w_triplets.begin(), w_triplets.end());
    }

    if(debug) {
        std::cout << "Weight diagonal matrix W (sparse):\n";
        std::cout << "W has " << W.nonZeros() << " nonzeros\n";

        // print the W matrix for debugging
        std::cout << "W matrix (first 5 rows and columns):\n";
        int rows = W.rows(), cols = W.cols();
        for (int i = 0; i < std::min(5, rows); ++i) {
            for (int j = 0; j < std::min(5, cols); ++j) {
                std::cout << W.coeff(i, j) << " ";
            }
            std::cout << "\n";
        }
    }
}

Eigen::MatrixXd ElectricalFlowNaive::buildPseudoInverse(Eigen::MatrixXd& L) {
    // 3. Compute pseudoinverse L†
    //    Note: L is symmetric, but singular (rank n-1)
    //    So we use Eigen's SelfAdjointEigenSolver
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(L);
    Eigen::VectorXd evals = eig.eigenvalues();
    Eigen::MatrixXd evecs = eig.eigenvectors();

    Eigen::MatrixXd Linv = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        if (evals[i] > 1e-8) { // avoid zero eigenvalue (nullspace of L)
            Linv += (1.0 / evals[i]) * (evecs.col(i) * evecs.col(i).transpose());
        }
    }


    if(debug) {
        std::cout << "Eigenvalues of L:\n" << evals.transpose() << "\n";
        std::cout << "Pseudo-inverse of L (dense):\n" << Linv << "\n";
    }

    return Linv;
}


Eigen::SparseMatrix<double> ElectricalFlowNaive::getFinalRoutingMatrix() {
    return M_avg; // Return the average M matrix as the final routing matrix
}


std::unordered_map<std::pair<int, int>, Eigen::VectorXd> ElectricalFlowNaive::getRoutingForCommodity(const std::vector<std::pair<int, int> >& _commodity) {
    std::unordered_map<std::pair<int, int>, Eigen::VectorXd>  routing;

    for(const auto& f : _commodity) {
        Eigen::VectorXd commodity_routing = Eigen::VectorXd::Zero(m_graph.getNumNodes());
        commodity_routing[f.first] = 1.0; // Set the flow for the first node
        commodity_routing[f.second] = -1.0; // Set the flow for the second node

        auto x = M_avg * commodity_routing; // Compute the routing for the commodity using the average M matrix

        routing[f] = x;
    }
    return routing; // Return the routing for the commodities
}


std::vector<FlowPath> ElectricalFlowNaive::decomposeFlowToPaths(
        int source,
        int sink,
        const Eigen::VectorXd& edge_flow
) {


    // Step 1: Build adjacency list with positive flow edges
    std::unordered_map<int, std::vector<std::pair<int, int>>> adj;
    std::unordered_map<std::pair<int, int>, double> flow_map;

    for (int e = 0; e < m; ++e) {
        auto [u, v] = edges[e];
        double f = edge_flow[e];
        if (f > 1e-8) {
            adj[u].emplace_back(v, e);
            flow_map[{u, v}] = f;
        }
    }

    std::vector<FlowPath> result;

    // Step 2: Decompose paths greedily
    while (true) {
        // BFS to find an s-t path
        std::unordered_map<int, std::pair<int, int>> parent; // v → (u, edge_id)
        std::queue<int> q;
        q.push(source);

        while (!q.empty() && !parent.count(sink)) {
            int u = q.front(); q.pop();
            for (auto [v, e_id] : adj[u]) {
                if (!parent.count(v)) {
                    parent[v] = {u, e_id};
                    q.push(v);
                }
            }
        }

        if (!parent.count(sink)) {
            break; // No more s-t paths
        }

        // Reconstruct path and find bottleneck
        std::vector<std::pair<int, int>> path;
        std::vector<int> edge_ids;
        double min_flow = std::numeric_limits<double>::infinity();
        int curr = sink;

        while (curr != source) {
            auto [prev, e_id] = parent[curr];
            path.emplace_back(prev, curr);
            edge_ids.push_back(e_id);
            min_flow = std::min(min_flow, flow_map[{prev, curr}]);
            curr = prev;
        }

        std::reverse(path.begin(), path.end());
        std::reverse(edge_ids.begin(), edge_ids.end());

        // Store this path and its flow
        result.push_back({path, min_flow});

        // Subtract flow along the path
        for (int e_id : edge_ids) {
            auto [u, v] = edges[e_id];
            flow_map[{u, v}] -= min_flow;
            if (flow_map[{u, v}] < 1e-8) {
                // Remove edge from residual graph
                auto& nbrs = adj[u];
                nbrs.erase(std::remove_if(nbrs.begin(), nbrs.end(),
                                          [&](const std::pair<int, int>& p) { return p.first == v; }), nbrs.end());
            }
        }
    }

    return result;
}


std::unordered_map<std::pair<int, int>, std::vector<FlowPath>>
ElectricalFlowNaive::getRoutingPathsForCommodity(const std::vector<std::pair<int, int>>& _commodity) {
    std::unordered_map<std::pair<int, int>, std::vector<FlowPath>> routing_paths;

    for (const auto& f : _commodity) {
        // Set up unit demand vector for this commodity
        Eigen::VectorXd commodity_routing = Eigen::VectorXd::Zero(m_graph.getNumNodes());
        commodity_routing[f.first] = 1.0;
        commodity_routing[f.second] = -1.0;

        // Compute electrical flow using M_avg
        Eigen::VectorXd x = M_avg * commodity_routing;

        // Print edge flow vector
        std::cout << "Raw edge flow for commodity (" << f.first << ", " << f.second << "):\n";
        std::cout << x.transpose() << "\n";

        // Run decomposition
        auto paths = decomposeFlowToPaths(f.first, f.second, x);

        // Print paths
        std::cout << "Decomposed flow paths for commodity (" << f.first << ", " << f.second << "):\n";
        for (const auto& p : paths) {
            std::cout << "  Flow " << p.amount << ": ";
            for (const auto& [u, v] : p.path) {
                std::cout << u << " → " << v << "  ";
            }
            std::cout << "\n";
        }

        routing_paths[f] = paths;
    }

    return routing_paths;
}


double ElectricalFlowNaive::getMaximumCongestion() {
    double max_cong = 0.0;

    std::vector<double> total_edge_congestion(edges.size(), 0);

    // also add the flow for the fixed vertex x_fixed
    for(int i = 0; i<n; i++) {
        for(int j = i+1; j<n; j++) {
            for(int edge_id = 0; edge_id < edges.size(); edge_id++) {

                double f_e_i = f_e_u.count({edge_id, i}) ? f_e_u[{edge_id, i}] : 0;
                double f_e_j = f_e_u.count({edge_id, j}) ? f_e_u[{edge_id, j}] : 0;

                double flow = f_e_i - f_e_j; // Get the flow for the edge
                if (std::abs(flow) > 1e-16) { // Only consider significant flows
                    f_e_st[{edge_id, i, j}] += flow; // Sum the flow for the commodity pair (i, j)
                }
            }
        }
    }

    if(debug) {
        // print out the flows
        for (auto &[edge_commo, flow]: f_e_st) {
            int edge_id = std::get<0>(edge_commo);
            int source = std::get<1>(edge_commo);
            int target = std::get<2>(edge_commo);
            std::cout << "Flow for edge (" << edges[edge_id].first << ", " << edges[edge_id].second
                      << ") commodity (" << source << " -> " << target << "): " << flow << "\n";
        }
    }


    // iterate over every commodity
    // and sum up the total flow per edge
    for(auto& [edge, flow] : f_e_st) {
        int edge_id = std::get<0>(edge);
        total_edge_congestion[edge_id] += std::abs(flow)/c_edges2capacities[edges[edge_id]]; // Sum up the absolute flow values for each edge
    }

    for(auto& congestion: total_edge_congestion) {
        if(congestion > max_cong) {
            max_cong = congestion; // Update the maximum congestion
        }
    }

    return max_cong;

}


double ElectricalFlowNaive::getCongestion(DemandMap& demands) {

    // TODO: remove the call getMaximumCongestion() here
    // as of now we only compute the congestion after the flow is computed
    getMaximumCongestion(); // Ensure the flow is computed before calculating congestion

    double max_cong = 0.0;

    std::vector<double> total_edge_congestion(edges.size(), 0);
    // iterate over every commodity
    for (auto &[edge_commo, flow]: f_e_st) {
        int edge_id = std::get<0>(edge_commo);
        int source = std::get<1>(edge_commo);
        int target = std::get<2>(edge_commo);

        for(auto& [d, demand_value] : demands) {
            if(source == d.first && target == d.second) {
                // If the edge is part of the current commodity
                total_edge_congestion[edge_id] += std::abs(flow)*demand_value/c_edges2capacities[edges[edge_id]]; // Sum up the absolute flow values for each edge
            }
        }
    }

    double max_cong_edge = *std::max_element(total_edge_congestion.begin(), total_edge_congestion.end());
    return max_cong_edge; // Return the maximum congestion across all edges
}