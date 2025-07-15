//
// Created by Mert Biyikli on 09.07.25.
//

#include "electrical_flow_naive.h"
#include <random>

void ElectricalFlowNaive::init(const RaeckeGraph &g) {
    m_graph = g;
    roh = std::sqrt(2.0*static_cast<double>(m_graph.getNumEdges())); // Initialize roh based on the number of nodes
    alpha_local = std::log2(m_graph.getNumNodes())*std::log2(m_graph.getNumNodes()); // Initialize alpha_local based on the number of nodes

    // Initialize the graph distances and capacities
    for (int v = 0; v < m_graph.getNumNodes(); ++v) {
        for (int u : m_graph.neighbors(v)) {
            x_edge2distance[{v, u}] = 1.0; // Initialize distances to 1.0
        }
    }
    this->cap_X = m_graph.getNumEdges();
    this->number_of_iterations = 8*roh*std::log(m_graph.getNumEdges())/alpha_local;

    extract_edge_list();
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

            double w = 1.0/(p + inv_m); // Calculate the weight based on the probability and inverse of the number of edges
            w_edges2weights[{v, u}] = w;
        }
    }
}
/*
Eigen::SparseMatrix<double> ElectricalFlowNaive::getRoutingMatrix() {
    // Build the incidence matrix B and the weight vector w
    Eigen::SparseMatrix<double> B = this->buildIncidence(m_graph.getNumNodes(), edges);
    LaplacianSolver solver;
    solver.init(m_graph, w_edges2weights, false);
    auto L = solver.GetLaplacianMatrixFromGraph(m_graph, w_edges2weights, false);

    // Build the diagonal weight matrix W
    Eigen::SparseMatrix<double> W = buildWeightDiag(weights);

    // Compute the routing matrix R = Bᵀ W B
    Eigen::SparseMatrix<double> R = B.transpose() * W * B;

    // compute the pseudoinsert of R


    return R;
}
*/

Eigen::SparseMatrix<double> ElectricalFlowNaive::buildWeightDiag(const std::vector<double>&) {
    int n = m_graph.getNumEdges();
    Eigen::SparseMatrix<double> W(n, n);

    return W;
}

void ElectricalFlowNaive::run() {
    inv_m = 1.0 / static_cast<double>(m_graph.getNumEdges());

    for(int t = 0; t < number_of_iterations; ++t) {


        updateEdgeDistances();

        // 2) Refresh our weights[] array in the same edge‐order
        for(int i = 0; i < cap_X; ++i) {
            auto e = edges[i];
            weights[i] = w_edges2weights[e];
        }

        //getRoutingMatrix();

        // 3) Build M = W B (Bᵀ W B)^{-1} via m solves
        Eigen::SparseMatrix<double> M;
        computeM_via_solves(M);

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

        auto aload = getApproxLoad();
        // print approxiload
        std::cout << "Approximate load:\n";
        for (const auto& kv : aload) {
            auto [edge, load] = kv;
            std::cout << "Edge (" << edge.first << ", " << edge.second << "): Load = " << load << "\n";
        }
        break;
    }
}


std::unordered_map<std::pair<int, int>, double> ElectricalFlowNaive::getApproxLoad() {

    std::unordered_map<std::pair<int, int>, double> approx_load;
    // build the incidence matrix B and the weight vector w
    Eigen::SparseMatrix<double> B = this->buildIncidence(m_graph.getNumNodes(), edges);
    LaplacianSolver solver;
    solver.init(m_graph, w_edges2weights, false);
    int m = m_graph.getNumEdges();

    // Get the Laplacian matrix L from the solver
    // todo
    auto C = getSketchMatrix(m_graph.getNumEdges(), m_graph.getNumNodes(), 0.5);

    auto X = B.transpose() * C.transpose();

    std::vector<Eigen::VectorXd> U_i;
    // for each column of X, solve the linear system with the Laplacian of the graph itself
    for(int i = 0; i<X.cols(); ++i) {
        Eigen::VectorXd x = X.col(i);
        Eigen::VectorXd u = solver.solve(x);
        U_i.push_back(u);
    }

    // Now we have U_i, which is a vector of vectors, where each vector corresponds to a column of X
    // cumulate the U_i vectors to get an matrix to transpose
    Eigen::MatrixXd U;
    for(int i = 0; i < U_i.size(); ++i) {
        U.row(i) = U_i[i];
    }

    // for each edge (u,v), compute the approximate load
    for(int i = 0; i < edges.size(); ++i) {
        auto b = B.row(i);

        double norm = recoverNorm(U.transpose(), b);

        approx_load[edges[i]] = w_edges2weights[edges[i]] * norm;

    }


    return approx_load;
    // get the approxi load from the recover norm
}

// Build the oriented incidence matrix B (E×V).  Each edge e=(u,v):
//   row e has +1 in column u, –1 in column v.
Eigen::SparseMatrix<double> ElectricalFlowNaive::buildIncidence(int V,std::vector<std::pair<int,int>>& _edges){
    int E = _edges.size();
    Eigen::SparseMatrix<double> B(E, V);
    std::vector<Eigen::Triplet<double>> T;
    T.reserve(2*E);
    for(int e = 0; e < E; ++e){
        auto [u,v] = _edges[e];
        T.emplace_back(e, u, +1.0);
        T.emplace_back(e, v, -1.0);
    }
    B.setFromTriplets(T.begin(), T.end());
    return B;
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
    assert(int(edges.size()) == m_graph.getNumEdges());
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
    std::cout << "Sketch matrix C (first 5 rows and columns):\n";
    int rows = C.rows(), cols = C.cols();
    for (int i = 0; i < std::min(5, rows); ++i) {
        for (int j = 0; j < std::min(5, cols); ++j) {
            std::cout << C(i, j) << " ";
        }
        std::cout << "\n";
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
