//
// Created by Mert Biyikli on 30.09.25.
//

#include "electrical_flow_optimized.h"
#include "../utils/my_math.h"
#include <random>
/*
    Main algorithm implementation
 */

void ElectricalFlowOptimized::init(bool debug, const std::string& solver_name)
{

    n = graph.getNumNodes();
    m = graph.getNumEdges()/2;


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
        double cap = graph.getEdgeCapacity(e);      // undirected capacity accessor
        edge_capacities[e]  = cap;
        edge_probabilities[e] = edge_distances[e] / cap_X;
        edge_weights[e] = std::pow(cap, 2) / (edge_probabilities[e] + inv_m);
    }


}


void ElectricalFlowOptimized::extract_edge_list() {

    // adjancency list version

    // -> for the adjacency version we may not need an extra vector edges anymore...
    edges.reserve(m);
    for (int e = 0; e < graph.getNumEdges(); e++) {
        auto [u, v] = graph.edgeEndpoints(e);
        if (u < v) {
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




void ElectricalFlowOptimized::run(LinearRoutingTable &table) {
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
            pure_oracles_running_times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(end_pure-start_pure).count());

            // accumulate edge flows for commodity (u -> x_fixed):
            // f_e(u) = w_e * (x[u]-x[v])

            // TODO: Make this iterating over all UNORDERED edges instead of ordered edges only
            for (int e = 0; e < m; ++e) {
                double fval = edge_weights[e] * (x[edges[e].first] - x[edges[e].second]);

                if (u == edges[e].first)      div_accum[u] += fval;
                else if (u == edges[e].second) div_accum[u] -= fval;
                if (std::abs(fval) > EPS) {
                    // if value is negative, we push along the anti-edge
                    const auto& [head, tail] = edges[e];
                    const int& original_edge_id = graph.getEdgeId(head, tail);


                    // std::cout << "electrical e: " << e << "\n";
                    if (fval < 0) {
                        int anti_edge = graph.getAntiEdge(original_edge_id);
                        if (anti_edge) {
                            // std::cout << "anti_e= " << anti_edge << " head=" << head << " tail=" << tail << "\n";
                            // std::cout << "found anti edge adding: " << std::abs(fval) << " for source " << u << "\n";
                        }
                        table.addFlow(anti_edge, u, std::abs(fval));
                    }else {
                        if (original_edge_id) {
                            // std::cout << "orig_e= " << original_edge_id << " head=" << head << " tail=" << tail << "\n";
                            // std::cout << "found anti edge adding: " << fval << " for source " << u << "\n";
                        }
                        table.addFlow(original_edge_id, u, std::abs(fval));
                    }
                }
            }

            for (int e = 0; e < m; ++e) {
                double fval = edge_weights[e] * (x[edges[e].second] - x[edges[e].first]);

                if (u == edges[e].second)      div_accum[u] += fval;
                else if (u == edges[e].first) div_accum[u] -= fval;
                if (std::abs(fval) > EPS) {

                    // if value is negative, we push along the anti-edge
                    auto [head, tail] = edges[e];
                    std::swap(head, tail);
                    const int& original_edge_id = graph.getEdgeId(head, tail);
                    // std::cout << "electrical e: " << e << "\n";
                    if (fval < 0) {
                        int anti_edge = graph.getAntiEdge(original_edge_id);
                        if (anti_edge) {
                            // std::cout << "anti_e= " << anti_edge << " head=" << head << " tail=" << tail << "\n";
                            // std::cout << "found anti edge adding: " << std::abs(fval) << " for source " << u << "\n";
                        }
                        table.addFlow(anti_edge, u, std::abs(fval));
                    }else {
                        if (original_edge_id) {

                            // std::cout << "orig_e= " << original_edge_id << " head=" << head << " tail=" << tail << "\n";
                            // std::cout << "found anti edge adding: " << fval << " for source " << u << "\n";
                        }
                        table.addFlow(original_edge_id, u, std::abs(fval));
                    }

                }
            }

            // for debugging purposes, print the flows
            for (int i = 0; i<table.src_ids.size(); i++) {
                for (int j = 0; j<table.src_ids[i].size(); j++) {
                    int source = table.src_ids[i][j];
                    double flow = table.src_flows[i][j];
                    if (flow != 0.0) {
                        auto [head, tail] = graph.edgeEndpoints(i);
                        // std::cout << "Edge " << head << " / " << tail << " has flow " << flow / iteration_count << " for source " << source << "\n";
                    }
                }
            }

        }


        getApproxLoad(load);
        updateEdgeDistances(load); // keeps 'weights' vector in sync

        auto end = std::chrono::high_resolution_clock::now();

        oracle_running_times.push_back(
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        );

    }


}


void ElectricalFlowOptimized::scaleFlowDown(LinearRoutingTable& table) {

    // scale the flow from the adjacency list flow
    if (iteration_count > 0) {
        const double inv_iters = 1.0 / static_cast<double>(iteration_count);
        for (int e = 0; e < graph.getNumEdges(); ++e) // dont use m here. m is undirected edges only
            for (double &val : table.src_flows[e]) val *= inv_iters;
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
        pure_oracles_running_times.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count());

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