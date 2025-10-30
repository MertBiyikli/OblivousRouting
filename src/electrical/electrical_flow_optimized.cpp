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
    this->debug = debug;

    // set algorithm parameters
    roh = std::sqrt(2.0*static_cast<double>(m)); // Initialize roh based on the number of nodes
    alpha_local = std::log2(n)*std::log2(n); // Initialize alpha_local based on the number of nodes
    this->cap_X = m;
    this->iteration_count = 8.0*roh*std::log(m)/alpha_local;
    this->inv_m = 1.0 / static_cast<double>(m);
    this->div_accum.assign(n, 0.0);


    // fix a node x
    this->x_fixed = 0; // Randomly select a fixed node x from the graph
    //tree_parent = buildBFSTree(x_fixed);


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
    // --- 1) stable MWU in log-space with damping & clipping ---
    const double eta   = 1.0 / (4.0 * roh);   // smaller than your 1/(2*roh)
    const double gclip = 10.0;                // gradient clip (tune 5–20)

    // Kahan sum for cap_X
    double capX = 0.0, c = 0.0;

    for (int e = 0; e < m; ++e) {
        double g = std::clamp(load[e], -gclip, gclip);
        // multiplicative: x_e <- x_e * exp(eta * g)
        double xe = edge_distances[e];
        if (xe > 0.0 && std::isfinite(g)) {
            double upd = std::exp(eta * g);
            xe *= upd;
            // avoid denormals
            if (!std::isfinite(xe)) xe = (g > 0 ? 1e308 : 1e-308);
            edge_distances[e] = xe;
        }
        // Kahan accumulate cap_X
        double y = xe - c;
        double t = capX + y;
        c = (t - capX) - y;
        capX = t;
    }
    cap_X = capX;

    // --- 2) recompute probabilities and weights with caps ---
    // floors / caps to keep W well-conditioned
    const double p_floor = 1.0 / (10.0 * m);       // soft floor
    const double w_floor = 1e-12;
    const double w_cap   = 1e12;                   // tune as needed

    for (int e = 0; e < m; ++e) {
        double x = edge_distances[e];
        if (x <= 0.0 || !std::isfinite(x)) x = 1e-12;

        double p = x / cap_X;
        if (!std::isfinite(p)) p = 1.0 / m;
        p = std::max(p, p_floor);                  // floor the prob

        edge_probabilities[e] = p;

        const double cap = edge_capacities[e];
        double w = (cap * cap) / (p + inv_m);
        if (!std::isfinite(w)) w = w_cap;
        w = std::min(std::max(w, w_floor), w_cap); // clamp
        edge_weights[e] = w;
    }

    // Tell AMG we're only changing numeric values (no rebuild here)
    amg->updateAllEdges(edge_weights, edges);  // ideally becomes value-only update
    amg->updateSolver();                       // later: make this infrequent (every ~25 iters)

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
            auto start_pure = std::chrono::high_resolution_clock::now();
            Eigen::VectorXd x = amg->solve(rhs);
            auto end_pure = std::chrono::high_resolution_clock::now();
            pure_oracle_running_times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(end_pure-start_pure).count());

            // accumulate edge flows for commodity (u -> x_fixed):
            // f_e(u) = w_e * (x[u]-x[v])
            for (int e = 0; e < m; ++e) {
                double fval = edge_weights[e] * (x[edges[e].first] - x[edges[e].second]);

                if (u == edges[e].first)      div_accum[u] += fval;
                else if (u == edges[e].second) div_accum[u] -= fval;
                if (std::abs(fval) > 1e-32) {
                    addFlow(e, u, fval); // per-adjacency list version
                    //f_e_u[{e, u}] += fval;
                }
            }
        }

        // (Optional) debug: print flows per (edge, source)
        if (debug) {
            printFlow();
        }

        getApproxLoad(load);
        updateEdgeDistances(load); // keeps 'weights' vector in sync

        auto end = std::chrono::high_resolution_clock::now();
        oracle_running_times.push_back(
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        );
    }


    scaleFlowDown();
    // enforceConservation();
    if (debug) printFlow();
}


void ElectricalFlowOptimized::scaleFlowDown() {

    // scale the flow from the adjacency list flow
    if (iteration_count > 0) {
        const double inv_iters = 1.0 / static_cast<double>(iteration_count);
        for (int e = 0; e < m; ++e)
            for (double &val : adj_f_e_u[e]) val *= inv_iters;
    }


    if ( debug ) {
        printFlow();
    }

    // do the recaling for the adjacency matrix as well
    std::vector<double> outgoingflow_f_s_x_fixed(n, 0);
    for (int e = 0; e < m; ++e) {
        const int u = edges[e].first;
        const int v = edges[e].second;
        const auto &ids  = adj_f_e_u_id[e];
        const auto &vals = adj_f_e_u[e];
        const int L = static_cast<int>(ids.size());

        for (int k = 0; k < L; ++k) {
            int s = ids[k];
            double f = vals[k];

            if ( f >= 0 ) {
                // direction u -> v
                if (s == u) {
                    outgoingflow_f_s_x_fixed[s] += f;
                }
            } else {
                // direction v -> u
                if (s == v) {
                    outgoingflow_f_s_x_fixed[s] += std::abs(f);
                }
            }
        }
    }

    // print outgoing flow
    for (int i = 0; i < n; ++i) {
        if (i == x_fixed) continue;
        double val = outgoingflow_f_s_x_fixed[i];
        if (debug ) std::cout << "outgoing flow: " << val << " for source " << i << "\n";
    }

    // Print out the div_accum values for debugging
    if (debug) {
        for ( int s = 0; s < n; ++s) {
            if (s == x_fixed) continue;
            std::cout << "div_accum[" << s << "] = " << div_accum[s] << "\n";
        }
    }

    // ---------- 3) Rescale all (e, s) by the commodity's net outflow ----------
    // After this, each commodity s should have unit net outflow at s.
    constexpr double EPS = 1e-12;
    for (int e = 0; e < m; ++e) {
        auto &ids  = adj_f_e_u_id[e];
        auto &vals = adj_f_e_u[e];
        const int L = static_cast<int>(ids.size());

        for (int k = 0; k < L; ++k) {
            const int s = ids[k];
            if (s == x_fixed) continue;  // no commodity starting at x_fixed

            const double denom = outgoingflow_f_s_x_fixed[s];
            if (std::abs(denom) <= EPS) {
                if (debug) {
                    std::cerr << "[scaleFlowDown] Warning: net_outflow ~ 0 for source "
                              << s << " (|denom|=" << std::abs(denom) << "), skipping rescale.\n";
                }
                continue;
            }

            const double scaled = vals[k] / denom;
            vals[k] = scaled;

            // Keep the old map in sync for downstream code (e.g., getFlowForCommodity).
            // The map key is (edge_id, source).
            //f_e_u[{e, s}] = scaled;
        }
    }
}

void ElectricalFlowOptimized::getApproxLoad(std::vector<double>& load) {
    const int ell = X.cols();
    // Preconditions: edge_wc precomputed, load sized to m, edge_* sized to m.

    // 1) Flat buffer for edge-major diffs
    edge_diffs.resize(static_cast<size_t>(m) * ell); // member or static thread-local if re-used

    Eigen::VectorXd rhs(n), sol(n);

    // 2) Fill diffs (AMG solve → edge diffs)
    for (int i = 0; i < ell; ++i) {
        rhs = X.col(i);

        auto start = std::chrono::high_resolution_clock::now();
        sol = amg->solve(rhs);  // prefer in-place if available
        auto end = std::chrono::high_resolution_clock::now();
        pure_oracle_running_times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count());

        const double* __restrict psol = sol.data();

        // stride on store (edge-major), but gives fast medians later
        for (int e = 0; e < m; ++e) {
            edge_diffs[static_cast<size_t>(e) * ell + i] =
                std::fabs(psol[edges[e].first] - psol[edges[e].second]);
        }
    }

    // 3) Medians + load
    #pragma omp parallel for schedule(static)
    for (int e = 0; e < m; ++e) {
        double* __restrict arr = &edge_diffs[static_cast<size_t>(e) * ell];
        std::nth_element(arr, arr + (ell >> 1), arr + ell);
        const double med = arr[ell >> 1];
        load[e] = edge_weights[e] * med;
    }

    // Optional: debug after, not inside hot loops
    if (debug) {
        for (int e = 0; e < m; ++e) {
            // print sampled edges only, or aggregate stats
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
    for (int e = 0; e<adj_f_e_u_id.size(); e++) {
        int u = edges[e].first;
        int v = edges[e].second;
        for (int s = 0; s<n; ++s) {
            for (int t = 0; t<n; t++) {
                if ( s == t ) continue;


                double flow(0);
                auto it = std::find(adj_f_e_u_id[e].begin(), adj_f_e_u_id[e].end(), s);
                if (it != adj_f_e_u_id[e].end()) {
                    int idx = std::distance(adj_f_e_u_id[e].begin(), it);
                    flow = adj_f_e_u[e][idx];
                }
                else {
                    auto it2 = std::find(adj_f_e_u_id[e].begin(), adj_f_e_u_id[e].end(), t);
                    if (it2 != adj_f_e_u_id[e].end()) {
                        int idx = std::distance(adj_f_e_u_id[e].begin(), it2);
                        flow = -adj_f_e_u[e][idx];
                    }
                }
                if (std::abs(flow) > 1e-9) {
                    // since this linear
                    if (flow <0) {
                        f_e_st[{v, u}][{s, t}] = -flow; // Store the reverse flow as well
                    }else {
                        f_e_st[{u, v}][{s, t}] = flow;
                    }
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

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(l) * _m);

    for (int i = 0; i < l; ++i) {
        for (int j = 0; j < _m; ++j) {
            triplets.emplace_back(i, j, dist(gen));
        }
    }

    Eigen::SparseMatrix<double> C(l, _m);
    C.setFromTriplets(triplets.begin(), triplets.end());
    C.makeCompressed();
    return C;
}