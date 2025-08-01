//
// Created by Mert Biyikli on 09.07.25.
//

#include "electrical_flow_naive.h"
#include <random>

void ElectricalFlowNaive::init(const RaeckeGraph &g, bool debug) {
    m_graph = g;

    n = m_graph.getNumNodes();
    m = m_graph.getNumEdges();

    // set algorithm parameters
    roh = std::sqrt(2.0*static_cast<double>(m_graph.getNumEdges())); // Initialize roh based on the number of nodes
    alpha_local = std::log2(m_graph.getNumNodes())*std::log2(m_graph.getNumNodes()); // Initialize alpha_local based on the number of nodes
    this->cap_X = m_graph.getNumEdges();
    this->number_of_iterations = 8.0*roh*std::log(m_graph.getNumEdges())/alpha_local;

    // initialize distance, weights, probabilities
    initEdgeDistances();
    extract_edge_list();

    // initialize the Laplacian Solver
    solver.init(m_graph, w_edges2weights); // Initialize the Laplacian solver with the current weights
    solver.setDebug(debug);
    this->debug = debug;
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
    for(int u = 0; u < m_graph.getNumNodes(); ++u) {
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
    inv_m = 1.0 / static_cast<double>(m_graph.getNumEdges());

    for(int t = 0; t < number_of_iterations; ++t) {

        // 3) Build M = W B (Bᵀ W B)^{-1} via m solves
        Eigen::SparseMatrix<double> M = getRoutingMatrix();

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

            x_edge2distance[edge] *= (1+(1/(2*roh))*load); // Update the distance with the approximate load
            cap_X += x_edge2distance[edge]; // Update the total capacity
        }

        updateEdgeDistances();

        //reak;
    }


    // take the average of the Ms matrices
    M_avg = Ms[0];
    for(int i = 1; i < Ms.size(); ++i) {
        M_avg += Ms[i];
    }
    M_avg /= static_cast<double>(Ms.size());


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
}


std::unordered_map<std::pair<int, int>, double> ElectricalFlowNaive::getApproxLoad() {

    std::unordered_map<std::pair<int, int>, double> approx_load;
    // build the incidence matrix B and the weight vector w
    if(!B.nonZeros()) {
        buildIncidence();
    }

    int m = m_graph.getNumEdges();
    int n = m_graph.getNumNodes();

    // Get the Laplacian matrix L from the solver
    // todo
    auto C = getSketchMatrix(m, n, 0.5);

    auto X = B.transpose() * C.transpose();

    Eigen::MatrixXd U(n, X.cols() );
    // for each column of X, solve the linear system with the Laplacian of the graph itself
    for(int i = 0; i<X.cols(); ++i) {
        Eigen::VectorXd x = X.col(i);
        Eigen::VectorXd u = solver.solve(x);
        U.col(i) = u;
    }

    // for each edge (u,v), compute the approximate load
    for(int i = 0; i < edges.size(); ++i) {
        auto b = B.row(i);

        double norm = recoverNorm(U.transpose(), b.transpose());

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




/*
 * This is the unoptimized version.
 * To reach full performance potential, replace the linear system
 * solver with the one frmo Spielamn to reach O(E log V) time.
 */
Eigen::MatrixXd ElectricalFlowNaive::getSketchMatrix(int m, int n, double epsilon) {
    // ℓ = O(ε⁻² log(1/δ))
    double delta = 1/log(n);
    int l = int(std::ceil( (2.0/(epsilon*epsilon)) * std::log(1.0/delta) ));
    Eigen::MatrixXd C(l, m);

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::cauchy_distribution<double> dist(0.0, 1.0);

    for(int i = 0; i < l; ++i)
        for(int j = 0; j < m; ++j)
            C(i,j) = dist(gen);

    // print sketch matrix for debugging
    if(debug) {
        std::cout << "Sketch matrix C (first 5 rows and columns):\n";
        int rows = C.rows(), cols = C.cols();
        for (int i = 0; i < std::min(5, rows); ++i) {
            for (int j = 0; j < std::min(5, cols); ++j) {
                std::cout << C(i, j) << " ";
            }
            std::cout << "\n";
        }
    }
    return C;
}

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
    for(int v = 0; v < m_graph.getNumNodes(); ++v) {
        for(int u : m_graph.neighbors(v)) {
            double distance = x_edge2distance[{v, u}];
            if(distance <= 0.0 || std::isnan(distance)) continue;

            // Update the probability based on the new distance
            double p = distance/cap_X;
            p_edge2probability[{v, u}] = p;

            double w = std::pow(c_edges2capacities[{v, u}], 2) * 1.0/(p + inv_m); // Calculate the weight based on the probability and inverse of the number of edges
            w_edges2weights[{v, u}] = w;

            solver.NonSyncUpdateEdgeWeights({v, u}, w); // Update the weights in the Laplacian solver
        }
    }
    solver.SyncEdgeWeights(); // Synchronize the weights in the Laplacian solver

    // Refresh our weights[] array in the same edge‐order
    refreshWeightMatrix();
    for(int i = 0; i < edges.size(); ++i) {
        auto e = edges[i];
        weights[i] = w_edges2weights[e];
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




Eigen::SparseMatrix<double>  ElectricalFlowNaive::getRoutingMatrix() {

    buildIncidence();
    buildWeightDiag();


    // 2. Compute Laplacian L = B * W * Bᵗ
    Eigen::SparseMatrix<double> WB = W*B;          // (n × m) * (m × m) = n × m
/*
    Eigen::MatrixXd M(m, n);

    for (int j = 0; j < n; ++j) {
        Eigen::VectorXd ej = Eigen::VectorXd::Zero(n);
        ej[j] = 1.0;

        // Solve using Julia-based Laplacian solver
        Eigen::VectorXd u = solver.solve(ej);

        // Multiply by WB → column j of M
        M.col(j) = WB * u;
    }

    // Convert to sparse and return
    return M.sparseView();*/

    if(debug) {
        std::cout << "Weighted incidence matrix WB (sparse):\n";
        std::cout << "WB has " << WB.nonZeros() << " nonzeros\n";
        // print the WB matrix for debugging
        std::cout << "WB matrix (first 5 rows and columns):\n";
        int rows = WB.rows(), cols = WB.cols();
        for (int i = 0; i < std::min(5, rows); ++i) {
            for (int j = 0; j < std::min(5, cols); ++j) {
                std::cout << WB.coeff(i, j) << " ";
            }
            std::cout << "\n";
        }
    }

    Eigen::MatrixXd L = (B.transpose()*WB).toDense(); // (n × m) * (m × n) = n × n


    if(debug) {
        std::cout << "Laplacian matrix L (dense):\n" << L << "\n";
    }

    auto Linv = buildPseudoInverse(L); // Compute the pseudo-inverse of the Laplacian matrix


    // 4. Compute M = W * Bᵗ * L†
    Eigen::MatrixXd M = WB*Linv;//Eigen::MatrixXd(m, n);

    if (debug) {
        std::cout << "Routing matrix M (dense):\n" << M << "\n";
    }

    return M.sparseView(); // Convert to sparse matrix for efficiency



    // (3) Allocate dense M matrix (edges x nodes)
    Eigen::MatrixXd M_dense(m, n);

    // compute the pseudo inverse of the Laplacian matrix


/*
    // TODO: remove this
    std::vector<double> abs_u;
    std::vector<double> abs_Bu;
    for (int j = 0; j < n; ++j) {

        // Unit vector e_j
        Eigen::VectorXd ej = Eigen::VectorXd::Zero(n);
        ej[j] = 1.0;

        // Solve L u = e_j
        Eigen::VectorXd u = solver.solve(ej);

        if(debug) {
            // print out the u vector for debugging
            std::cout << "u vector for column " << j << ":\n";
            for (int i = 0; i < u.size(); ++i) {
                std::cout << "Node " << i << ": " << u[i] << "\n";
            }

            for (int i = 0; i < u.size(); ++i) {
                abs_u.push_back( u[i] );
            }
        }

        // Compute Bu
        Eigen::VectorXd Bu(m);
        for (int e = 0; e < m; ++e) {
            auto [i, k] = edges[e];
            Bu[e] = u[i] - u[k];
        }

        // Multiply by weights
        for (int e = 0; e < m; ++e) {
            Bu[e] *= W.coeff(e, e);
        }

        if(debug) {

            for (auto& i : Bu) {
                abs_Bu.push_back(i);
            }

            // print out the Bu vector for debugging
            std::cout << "Bu vector for column " << j << ":\n";
            for (int e = 0; e < Bu.size(); ++e) {
                std::cout << "Edge (" << edges[e].first << ", " << edges[e].second << "): " << Bu[e] << "\n";
            }

        }

        // Store as column j
        M_dense.col(j) = Bu;
    }

    // Convert dense to sparse if desired
    Eigen::SparseMatrix<double> M_sparse = M_dense.sparseView();

    return M_sparse;
*/
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
        std::cout << "Eigenvectors of L (first 5 columns):\n" << evecs.leftCols(5) << "\n";
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


double ElectricalFlowNaive::getMaximumCongestion() const {
    double max_cong = 0.0;

    std::vector<double> total_edge_flow(edges.size(), 0);
    // iterate over every commodity
    // and sum up the total flow per edge
    Eigen::VectorXd demands(n);
    for(int s = 0; s<n; s++) {
        for(int t = s+1; t<n; t++) {
            demands.setZero();
            demands[s] = 1.0;
            demands[t] = -1.0;

            auto flow = M_avg*demands;

            for(int e = 0; e<flow.size(); ++e) {
                total_edge_flow[e] += std::abs(flow[e]);
            }
        }
    }


    for(int edge_id = 0; edge_id < total_edge_flow.size(); edge_id++){
        double flow = total_edge_flow[edge_id];
        int u = edges[edge_id].first, v = edges[edge_id].second;

        double cap = ((u < v) ? m_graph.getEdgeCapacity(u, v) : m_graph.getEdgeCapacity(v, u));

        double cong_e = std::abs(flow)/cap;
        if(cong_e > max_cong) {
            max_cong = cong_e;
        }
    }
    return max_cong;
}