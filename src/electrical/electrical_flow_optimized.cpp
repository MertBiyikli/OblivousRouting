//
// Created by Mert Biyikli on 30.09.25.
//

#include "electrical_flow_optimized.h"
#include "../utils/my_math.h"
#include <random>
/*
    Main algorithm implementation
 */

void ElectricalFlowOptimized::init(const Graph &g, bool debug, const std::string& solver_name)
{
    m_graph = g;

    n = m_graph.getNumNodes();
    m = m_graph.getNumEdges();

    adj_f_e_u_id.resize(m);
    adj_f_e_u.resize(m);

    // set algorithm parameters
    roh = std::sqrt(2.0*static_cast<double>(m)); // Initialize roh based on the number of nodes
    alpha_local = std::log2(n)*std::log2(n); // Initialize alpha_local based on the number of nodes
    this->cap_X = m;
    this->iteration_count = 8.0*roh*std::log(m)/alpha_local;
    this->inv_m = 1.0 / static_cast<double>(m);


    // fix a node x
    this->x_fixed = rand()%n; // Randomly select a fixed node x from the graph

    initEdgeDistances();
    buildWeightDiag();
    buildIncidence();

    // init AMG
    amg = SolverFactory::create(solver_name);
    amg->init(edge_weights, n, edges, debug);

    // U: diagonal of capacities (per-edge)
    U = Eigen::SparseMatrix<double>(m, m);
    for (int e = 0; e < m; ++e) U.insert(e, e) = edge_capacities[e];

    SketchMatrix_T = getSketchMatrix(m, n, 0.5).transpose();

    Eigen::MatrixXd UCt = SketchMatrix_T; // m × ℓ
    for (int e = 0; e < m; ++e)
        UCt.row(e) *= U.coeff(e, e);

    auto X_dense = (B.transpose() * UCt); // n × ℓ
    X = X_dense.sparseView();

}

void ElectricalFlowOptimized::initEdgeDistances() {
    extract_edge_list();

    // note that the edges are stored undirected
    // --- Per-edge init ---
    edge_distances.assign(m, 1.0);
    edge_capacities.resize(m);
    edge_probabilities.resize(m);
    edge_weights.resize(m);

    for (int e = 0; e < m; ++e) {
        auto [u,v] = edges[e];
        double cap = m_graph.getEdgeCapacity(u, v);      // undirected capacity accessor
        edge_capacities[e]  = cap;
        edge_probabilities[e] = edge_distances[e] / cap_X;
        edge_weights[e] = std::pow(cap, 2) / (edge_probabilities[e] + inv_m);
    }


}


void ElectricalFlowOptimized::extract_edge_list() {

    // adjancency list version

    // -> for the adjacency version we may not need an extra vector edges anymore...
    edges.reserve(m);
    for (int u = 0; u<n; u++) {
        for (int idx = 0; idx < m_graph.neighbors(u).size(); ++idx) {
            const int v = m_graph.neighbors(u)[idx];
            if (u < v) {
                edges.emplace_back(u,v);
            }
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
    if (B.nonZeros() == 0) {
        B = Eigen::SparseMatrix<double>(m, n);
        std::vector<Eigen::Triplet<double>> T; T.reserve(2*m);
        for (int e = 0; e < m; ++e) {
            auto [u,v] = edges[e]; // u < v
            T.emplace_back(e, u, -1.0);
            T.emplace_back(e, v, +1.0);
        }
        B.setFromTriplets(T.begin(), T.end());
    }
}

void ElectricalFlowOptimized::buildWeightDiag() {
    if (!W.nonZeros()) {
        W = Eigen::SparseMatrix<double>(m, m);
        std::vector<Eigen::Triplet<double>> w_triplets; w_triplets.reserve(m);
        for (int e = 0; e < m; ++e)
            w_triplets.emplace_back(e, e, edge_weights[e]);
        W.setFromTriplets(w_triplets.begin(), w_triplets.end());
    }
}

void ElectricalFlowOptimized::refreshWeightMatrix() {
    if (!W.nonZeros()) {
        buildWeightDiag();
    } else {
        for (int e = 0; e < m; ++e)
            W.coeffRef(e, e) = edge_weights[e];
    }
}

void ElectricalFlowOptimized::updateEdgeDistances(const std::vector<double>& load) {

    // Update distances per-edge, then mirror to per-adj
    cap_X = 0.0;
    for (int e = 0; e < m; ++e) {
        if (load[e] <= 0.0 || std::isnan(load[e])) continue;
        edge_distances[e] *= (1 + (1/(2*roh)) * load[e]);
        cap_X += edge_distances[e];
    }

    // Per-edge recompute
    for (int e = 0; e < m; ++e) {
        double x = edge_distances[e];
        if (x <= 0.0 || std::isnan(x)) continue;

        double p = x / cap_X;
        edge_probabilities[e] = p;

        double cap = edge_capacities[e];
        double w   = std::pow(cap, 2) / (p + inv_m);
        edge_weights[e] = w;
    }
    amg->updateAllEdges(edge_weights, edges);
    amg->updateSolver();

    refreshWeightMatrix();
}



void ElectricalFlowOptimized::run() {
    // Sparse RHS with only two nonzeros per demand (+1 at u, -1 at x_fixed)
    Eigen::SparseVector<double> rhs(n);
    rhs.reserve(2);
    std::vector<double> load(m, 0.0);

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
            // f_e(u) = w_e * (x[u]-x[v])
            for (int e = 0; e < m; ++e) {
                double fval = edge_weights[e] * (x[edges[e].first] - x[edges[e].second]);
                if (std::abs(fval) > 1e-16) {
                    addFlow(e, u, fval); // per-adjacency list version
                    //f_e_u[{e, u}] += fval;
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

        getApproxLoad(load);
        updateEdgeDistances(load); // keeps 'weights' vector in sync

        auto end = std::chrono::high_resolution_clock::now();
        oracle_running_times.push_back(
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        );
    }
/*
    // check if the flow map is in sync with the adjacency list version
    bool in_sync = true;
    for (int e = 0; e < m; ++e) {
        for (int idx = 0; idx < adj_f_e_u_id[e].size(); ++idx) {
            int u = adj_f_e_u_id[e][idx];
            double f_adj = adj_f_e_u[e][idx];
            double f_map = f_e_u[{e, u}];
            if (std::abs(f_adj - f_map) > 1e-9) {
                in_sync = false;
                std::cerr << "Mismatch for edge_id=" << e << ", u=" << u
                          << ": adj_list=" << f_adj << ", map=" << f_map << "\n";
            }
        }
    }

    // check the reverse direction
    for(const auto& [e_u, f_map] : f_e_u) {
        int e = e_u.first;
        int u = e_u.second;
        auto it = std::find(adj_f_e_u_id[e].begin(), adj_f_e_u_id[e].end(), u);
        if (it == adj_f_e_u_id[e].end()) {
            in_sync = false;
            std::cerr << "Missing entry in adjacency list for edge_id=" << e << ", u=" << u
                      << ": map=" << f_map << "\n";
        } else {
            int idx = std::distance(adj_f_e_u_id[e].begin(), it);
            double f_adj = adj_f_e_u[e][idx];
            if (std::abs(f_adj - f_map) > 1e-9) {
                in_sync = false;
                std::cerr << "Mismatch for edge_id=" << e << ", u=" << u
                          << ": adj_list=" << f_adj << ", map=" << f_map << "\n";
            }
        }
    }*/
/*
    // Average flows across iterations
    for (auto &kv : f_e_u) {
        kv.second /= static_cast<double>(iteration_count);
    }*/

    //scaleFlowDown();
}


void ElectricalFlowOptimized::scaleFlowDown() {

    // scale the flow from the adjacency list flow
    for(int e = 0; e < m; ++e) {
        for(int idx = 0; idx < adj_f_e_u_id[e].size(); ++idx) {
            int u = adj_f_e_u_id[e][idx];
            double f_adj = adj_f_e_u[e][idx];
            f_adj /= static_cast<double>(iteration_count);
            adj_f_e_u[e][idx] = f_adj;
            f_e_u[{e, u}] = f_adj; // keep the map in sync
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
            if (debug)
                std::cerr << "Warning: No outgoing flow for source " << source << ". Skipping scaling.\n";
        }
    }
}

void ElectricalFlowOptimized::getApproxLoad(std::vector<double>& load) {
    int ell = X.cols();
    std::vector<std::vector<double>> edge_diffs(m);
    for (int e = 0; e < m; ++e) edge_diffs[e].reserve(ell);

    Eigen::VectorXd rhs(n), sol(n);
    for (int i = 0; i < ell; ++i) {
        rhs = X.col(i);
        sol = amg->solve(rhs);

        for (int e = 0; e < m; ++e) {
            auto [u,v] = edges[e];
            edge_diffs[e].push_back(std::abs(sol[u] - sol[v]));
        }
    }

    for (int e = 0; e < m; ++e) {
        auto &arr = edge_diffs[e];
        size_t mid = arr.size() / 2;
        std::nth_element(arr.begin(), arr.begin()+mid, arr.end());
        double med = arr[mid];
        load[e] = edge_weights[e] * med / edge_capacities[e];
        if (debug) {
            auto [u,v] = edges[e];
            std::cout << "Edge ("<<u<<","<<v<<") load="<<load[e]<<"\n";
        }
    }
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
// adjacency list version
    // TODO: uncomment this

    int e = 0;
    for (int u = 0; u<n; ++u) {
        for (int idx = 0; idx < m_graph.neighbors(u).size(); ++idx) {
            int v = m_graph.neighbors(u)[idx];
            if ( u > v ) continue;

            for (int s = 0; s<n; ++s) {
                for (int t = 0; t<n; t++) {
                    if ( s == t ) continue;

                    double flow = getFlowForCommodity(e, s, t);
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
            e++;
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