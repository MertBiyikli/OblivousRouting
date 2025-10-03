//
// Created by Mert Biyikli on 30.09.25.
//

#include "electrical_flow_optimized.h"
#include "../utils/my_math.h"
#include <random>
/*
    Main algorithm implementation
 */

void ElectricalFlowOptimized::init(const Graph &g,const std::string& solver_name,  bool debug)
{
    m_graph = g;

    n = m_graph.getNumNodes();
    m = m_graph.getNumEdges();

    // set algorithm parameters
    roh = std::sqrt(2.0*static_cast<double>(m)); // Initialize roh based on the number of nodes
    alpha_local = std::log2(n)*std::log2(n); // Initialize alpha_local based on the number of nodes
    this->cap_X = m;
    this->iteration_count = 8.0*roh*std::log(m)/alpha_local;
    this->inv_m = 1.0 / static_cast<double>(m);


    // fix a node x
    this->x_fixed = rand()%n; // Randomly select a fixed node x from the graph

    // initialize distance, weights, probabilities
    initEdgeDistances();
    extract_edge_list();

    // initialize the Laplacian Solver
    // Choose solver type dynamically (string or enum)
    amg = SolverFactory::create(solver_name);   // or "amg_bicgstab" / "eigen_direct"
    amg->init(w_edges2weights, n, debug);
    //amg->init(w_edge_weight, n,g, debug);
    this->debug = debug;


    // ToDo: remove this later on...


    // amg->check_openmp_runtime();
    // compute diagonal matrix consisting of only the capacities
    U = Eigen::SparseMatrix<double>(m, m);
    for(int i = 0; i < m; ++i) {
        U.insert(i, i) = c_edges2capacities[edges[i]]; // Insert the capacities into the diagonal matrix
    }

    this->SketchMatrix_T = getSketchMatrix(m, n, 0.5).transpose(); // Initialize the sketch matrix with epsilon = 0.5


    // Precompute X = Bᵀ U Cᵀ once
    buildIncidence();
    Eigen::MatrixXd UCt = SketchMatrix_T; // m × ℓ
    for (int e = 0; e < m; ++e)
        UCt.row(e) *= U.coeff(e,e);
    auto X_dense = (B.transpose() * UCt ); // n × ℓ
    X = X_dense.sparseView();

}

void ElectricalFlowOptimized::initEdgeDistances() {

/*
    this->x_edge_distance.resize(n);
    this->p_edge_probability.resize(n);
    this->w_edge_weight.resize(n);
    this->c_edge_capacity.resize(n);

    // Initialize the edge distances and probabilities
    for(int u = 0; u < n; ++u) {
        for(int v : m_graph.neighbors(u)) {
            double distance = 1.0;

            x_edge_distance[u].push_back(distance);
            double p_ = distance/cap_X;
            p_edge_probability[u].push_back(p_);
            double cap = m_graph.getEdgeCapacity(u, v);
            c_edge_capacity[u].push_back(cap);
            double w_ = std::pow(cap, 2) * 1.0/(p_ + inv_m);
            w_edge_weight[u].push_back(w_);

        }
    }
*/
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


void ElectricalFlowOptimized::extract_edge_list() {
    // we store only the undirected edges
    edges.reserve(m);
    for(auto &kv : x_edge2distance) {
        auto [u,v] = kv.first;
        if(u < v) {
            edges.emplace_back(u,v);
        }
    }



    // sort the edges based on the first node, then second node
    std::sort(edges.begin(), edges.end(),
              [](const std::pair<int,int> &a, const std::pair<int,int> &b) {
                  if (a.first != b.first) return a.first < b.first;
                  return a.second < b.second;
              });

}



void ElectricalFlowOptimized::buildIncidence()
{
    if(B.nonZeros() == 0) {
        B = Eigen::SparseMatrix<double>(m, n);
        std::vector<Eigen::Triplet<double>> T;
        T.reserve(2*m);
        for(int e = 0; e < m; ++e){
            // auto [u,v] = edges[e];
            T.emplace_back(e, edges[e].first, -1.0);
            T.emplace_back(e, edges[e].second, +1.0);
        }
        B.setFromTriplets(std::make_move_iterator(T.begin()),
                  std::make_move_iterator(T.end()));
    }
}

void ElectricalFlowOptimized::buildWeightDiag() {
    if(!W.nonZeros()) {
        // 1. Build W as a sparse diagonal matrix
        W = Eigen::SparseMatrix<double>(m, m);
        std::vector<Eigen::Triplet<double>> w_triplets(m);
        for (int e = 0; e < m; ++e) {
            double weight = w_edges2weights[edges[e]];
            w_triplets.emplace_back(e, e, weight); // Add a triplet for the diagonal entry
        }
        W.setFromTriplets(w_triplets.begin(), w_triplets.end());
    }

/*
    // TODO: NOT Yet finished -> unordered_map -> adj_list of edge weights
    if (!W.nonZeros()) {
        // 1. Build W as a sparse diagonal matrix
        W = Eigen::SparseMatrix<double>(m, m);
        std::vector<Eigen::Triplet<double>> w_triplets(m);
        int e_id = 0;
        for (int v = 0; v < n; v++) {
            for (int u_iter = 0; u_iter< m_graph.neighbors(v).size(); u_iter++) {
                if (v > m_graph.neighbors(v)[u_iter]) continue; // only add one direction of the edge
                // auto edge = std::make_pair(v, m_graph.neighbors(v)[u_iter]);

                // int e = std::distance(edges.begin(), std::find(edges.begin(), edges.end(), edge));
                double weight = w_edge_weight[v][u_iter];
                w_triplets.emplace_back(e_id, e_id, weight); // Add a triplet for the diagonal entry
                e_id++;
            }
        }
    }*/
}

void ElectricalFlowOptimized::updateEdgeDistances() {
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

        //amg->updateEdge(v, u, w);
        //solver.NonSyncUpdateEdgeWeights({v, u}, w); // Update the weights in the Laplacian solver

    }
    //solver.SyncEdgeWeights(); // Synchronize the weights in the Laplacian solver

    // TODO: change this by only changing the numeric values of the matrix solver
    //amg->init(w_edges2weights, n);
    amg->updateAllEdges(w_edges2weights);
    amg->updateSolver();

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

    // for now comment the below code out
    /*
    /// TODO: New:

    // TODO: here comes the implementation for using adj_list as edge weights
    cap_X = 0.0;

    for (int v = 0; v < n; v++) {
        for (int u_iter = 0; u_iter < m_graph.neighbors(v).size(); u_iter++) {
            int u = m_graph.neighbors(v)[u_iter];
            if (v > u) continue;

            double dist = x_edge_distance[v][u_iter];
            if (dist <= 0.0 || std::isnan(dist)) continue;

            double p_ = dist / cap_X;
            p_edge_probability[v][u_iter] = p_;

            double cap = m_graph.getEdgeCapacity(v, u);
            double w_ = std::pow(cap, 2) * 1.0/(p_ + inv_m);
            w_edge_weight[v][u_iter] = w_;

            // look up edge id and update edge vector
            int eid = edge_id_map[make_key(v,u)];
            edges[eid].weight = w_;
        }
    }


    refreshWeightMatrix();
    amg->init(w_edge_weight, n, m_graph);
    */
}

void ElectricalFlowOptimized::refreshWeightMatrix() {
    if (!W.nonZeros()) {
        W = Eigen::SparseMatrix<double>(m, m);
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(m);
        for (int e = 0; e < m; e++) {
            double w = w_edges2weights[edges[e]];
            triplets.emplace_back(e, e, w);
        }
        W.setFromTriplets(triplets.begin(), triplets.end());
    } else {
        for (int e = 0; e < m; e++) {
            double w = w_edges2weights[edges[e]];
            W.coeffRef(e,e) = w;
        }
    }
}


void ElectricalFlowOptimized::run() {
    // Sparse RHS with only two nonzeros per demand (+1 at u, -1 at x_fixed)
    Eigen::SparseVector<double> rhs(n);
    rhs.reserve(2);

    for (int t = 0; t < this->iteration_count; ++t) {
        auto start = std::chrono::high_resolution_clock::now();

        // --- main loop over sources (u -> x_fixed) ---
        for (int u = 0; u < n; ++u) {
            if (u == x_fixed) continue;

            // build RHS (sum is zero automatically)
            rhs.setZero();
            rhs.coeffRef(u)       =  1.0;
            rhs.coeffRef(x_fixed) = -1.0;

            // Solve L x = b (solver returns dense potentials; that’s expected)
            // If your amg wrapper doesn't accept SparseVector, use: amg->solve(rhs.toDense())
            Eigen::VectorXd x = amg->solve(rhs);

            // accumulate edge flows for commodity (u -> x_fixed):
            // f_e = w_e * (x[a] - x[b])
            // Use the contiguous 'weights' vector to avoid unordered_map in hot loop.
            for (int e = 0; e < m; ++e) {
                const auto [a, b] = edges[e];
                const double fval = w_edges2weights[edges[e]] * (x[a] - x[b]);
                if (std::abs(fval) > 1e-16) {
                    f_e_u[{e, u}] += fval;
                }
            }
        }

        // (Optional) debug: print flows per (edge, source)
        if (debug) {
            for (auto &kv : f_e_u) {
                const auto &edge_id_and_src = kv.first;
                const int edge_id = edge_id_and_src.first;
                const int src     = edge_id_and_src.second;
                const auto [a, b] = edges[edge_id];
                std::cout << "Flow e=(" << a << "," << b << ") commodity (" << src
                          << " -> " << x_fixed << "): " << kv.second << "\n";
            }
        }

        // Approximate load & weight update (unchanged)
        auto aload = getApproxLoad();

        // --- Update distances based on loads ---
        cap_X = 0.0;
        for (int e = 0; e < m; ++e) {
            auto edge = edges[e];
            double load = aload[e];
            if (load <= 0.0 || std::isnan(load)) continue;

            x_edge2distance[edge] *= (1 + (1/(2*roh)) * load);
            cap_X += x_edge2distance[edge];
        }


        updateEdgeDistances(); // keeps 'weights' vector in sync

        auto end = std::chrono::high_resolution_clock::now();
        oracle_running_times.push_back(
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        );
    }

    // Average flows across iterations
    for (auto &kv : f_e_u) {
        kv.second /= static_cast<double>(iteration_count);
    }

    scaleFlowDown();
}


void ElectricalFlowOptimized::scaleFlowDown() {
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
            if (debug)
                std::cerr << "Warning: No outgoing flow for source " << source << ". Skipping scaling.\n";
        }
    }
}

std::vector<double> ElectricalFlowOptimized::getApproxLoad() {
    int ell = X.cols(); // sketch dimension

    // Per-edge scratch storage
    std::vector<std::vector<double>> edge_diffs(m);
    for (int e = 0; e < m; ++e) edge_diffs[e].reserve(ell);

    Eigen::VectorXd rhs(n), sol(n);

    // Solve for each RHS column
    for (int i = 0; i < ell; ++i) {
        rhs = X.col(i);
        sol = amg->solve(rhs);   // ⬅️ reuse AMG hierarchy, only numeric update between runs

        // Accumulate differences for all edges
        // #pragma omp parallel for if(m > 200) schedule(static)
        for (int e = 0; e < m; ++e) {
            auto [u,v] = edges[e];
            double diff = std::abs(sol[u] - sol[v]);
            edge_diffs[e].push_back(diff);
        }
    }

    // Compute median per edge + weight/capacity scaling
    std::vector<double> approx_load(m);
// #pragma omp parallel for if(m > 200) schedule(static)
    for (int e = 0; e < m; ++e) {
        auto &arr = edge_diffs[e];
        size_t mid = arr.size()/2;
        std::nth_element(arr.begin(), arr.begin()+mid, arr.end());
        double med = arr[mid];
        approx_load[e] = w_edges2weights[edges[e]] * med / c_edges2capacities[edges[e]];

        if(debug) {
            const auto& [u,v] = edges[e];
            std::cout << "Edge ("<<u<<","<<v<<") load="<<approx_load[e]<<"\n";
        }
    }

    return approx_load;
}

// Compute the median absolute difference between two node potentials
// across all ℓ sampled solutions.
//   diffs_u = potentials[u][*]  (length ℓ)
//   diffs_v = potentials[v][*]  (length ℓ)
double ElectricalFlowOptimized::recoverNorm(const std::vector<double>& diffs_u,
                                        const std::vector<double>& diffs_v) {
    assert(diffs_u.size() == diffs_v.size());
    size_t L = diffs_u.size();

    std::vector<double> abs_diffs;
    abs_diffs.reserve(L);
    for (size_t i = 0; i < L; ++i) {
        abs_diffs.push_back(std::abs(diffs_u[i] - diffs_v[i]));
    }

    size_t mid = L / 2;
    std::nth_element(abs_diffs.begin(), abs_diffs.begin() + mid, abs_diffs.end());
    return abs_diffs[mid];
}


void ElectricalFlowOptimized::storeFlow(){
    for (int edge_id = 0; edge_id < edges.size(); edge_id++) {
        int u = edges[edge_id].first, v = edges[edge_id].second;
        for (int s = 0; s<n; ++s) {
            for (int t = 0; t<n; t++) {
                if ( s == t ) continue;

                double flow = getFlowForCommodity(edge_id, s, t);
                if (std::abs(flow) > 1e-9) {

                    // since this linear
                    if (flow <0) {
                        f_e_st[{v, u}][{s, t}] = -flow; // Store the reverse flow as well
                    }else {
                        f_e_st[{u, v}][{s, t}] = flow;
                    }
                    //f_e_st[{v, u}][{t, s}] = -flow; // Store the reverse flow as well
                }
            }
        }
    }
}

double ElectricalFlowOptimized::getFlowForCommodity(int edge_id, int source, int target) {
    if(target == x_fixed) {
        return - (f_e_u.count({edge_id, source}) ? f_e_u.at({edge_id, source}) : 0.0);
    }else if(source == x_fixed) {
        return f_e_u.count({edge_id, target}) ? f_e_u.at({edge_id, target}) : 0.0;
    }else{
        // Not a fixed-x commodity; you'll reconstruct from two x-based flows
        double f_s = f_e_u.count({edge_id, source}) ? f_e_u.at({edge_id, source}) : 0.0;
        double f_t = f_e_u.count({edge_id, target}) ? f_e_u.at({edge_id, target}) : 0.0;
        return f_t - f_s;
    }
}

Eigen::SparseMatrix<double> ElectricalFlowOptimized::getSketchMatrix(int _m, int _n, double eps) {
    double c = 1.1;
    double delta = NegativeExponent(_n, 10); // Set delta to a small value or a default value
    double epsilon = eps; // Set epsilon to the provided value or a default value

    int l = (c/(epsilon*epsilon) * std::log(1.0/delta) );

    if(debug) {
        std::cout << "Generating sketch matrix with l = " << l << ", m = " << _m << ", n = " << _n << ", epsilon = " << epsilon << "\n";
    }

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::cauchy_distribution<double> dist(0.0, 1.0);

    Eigen::SparseMatrix<double> C(l, _m);
    C.reserve(l * _m); // each row has _m nonzeros

    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < _m; ++j) {
            double rand_value = dist(gen);
            C.insert(i, j) = rand_value;
        }
    }
    C.makeCompressed();
    return C;
}