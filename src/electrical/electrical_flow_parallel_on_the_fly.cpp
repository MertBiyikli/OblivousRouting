//
// Created by Mert Biyikli on 22.10.25.
//

#include "electrical_flow_parallel_on_the_fly.h"

#include "../utils/my_math.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*
    Main algorithm implementation
 */

void ElectricalFlowParallelOnTheFly::init(const GraphADJ &g, bool debug, const std::string& solver_name)
{
    // print how many threads are running
#ifdef  _OPENMP
    std::cout << "Number of threads: " << omp_get_max_threads() << "\n";
#endif

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


    initEdgeDistances();
    buildWeightDiag();
    buildIncidence();

    // init AMG
    amg = SolverFactory::create("amg_parallel");
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

void ElectricalFlowParallelOnTheFly::initEdgeDistances() {
    extract_edge_list();

    // note that the edges are stored undirected
    // --- Per-edge init ---
    edge_distances.assign(m, 1.0);
    edge_capacities.resize(m);
    edge_probabilities.resize(m);
    edge_weights.resize(m);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int e = 0; e < m; ++e) {
        auto [u,v] = edges[e];
        double cap = m_graph.getEdgeCapacity(u, v);      // undirected capacity accessor
        edge_capacities[e]  = cap;
        edge_probabilities[e] = edge_distances[e] / cap_X;
        edge_weights[e] = std::pow(cap, 2) / (edge_probabilities[e] + inv_m);
    }


}


void ElectricalFlowParallelOnTheFly::extract_edge_list() {

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



void ElectricalFlowParallelOnTheFly::buildIncidence()
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

void ElectricalFlowParallelOnTheFly::buildWeightDiag() {
    if (!W.nonZeros()) {
        W = Eigen::SparseMatrix<double>(m, m);
        std::vector<Eigen::Triplet<double>> w_triplets; w_triplets.reserve(m);
        for (int e = 0; e < m; ++e)
            w_triplets.emplace_back(e, e, edge_weights[e]);
        W.setFromTriplets(w_triplets.begin(), w_triplets.end());
    }
}

void ElectricalFlowParallelOnTheFly::refreshWeightMatrix() {
    if (!W.nonZeros()) {
        buildWeightDiag();
    } else {
        for (int e = 0; e < m; ++e)
            W.coeffRef(e, e) = edge_weights[e];
    }
}

void ElectricalFlowParallelOnTheFly::updateEdgeDistances(const std::vector<double>& load) {
    // --- 1) stable MWU in log-space with damping & clipping ---

    // Update distances per-edge, then mirror to per-adj
    cap_X = 0.0;
#ifdef _OPENMP
    #pragma omp parallel for reduction(+:cap_X) schedule(dynamic)
#endif
    for (int e = 0; e < m; ++e) {
        if (load[e] <= 0.0 || std::isnan(load[e])) continue;
        edge_distances[e] *= (1 + (1/(2*roh)) * load[e]);
        cap_X += edge_distances[e];
    }

    // Per-edge recompute
#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
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




void ElectricalFlowParallelOnTheFly::run() {
    // Sparse RHS with only two nonzeros per demand (+1 at u, -1 at x_fixed)
    Eigen::SparseVector<double> rhs(n);
    rhs.reserve(2);
    std::vector<double> load(m, 0.0);

#ifdef _OPENMP
    int num_threads = omp_get_max_threads();
#else
    int num_threads = 1;
#endif

    const int batch_size = m/num_threads;   // tune: 4–16 works best
    const int chunk       = batch_size;

    std::vector<std::vector<std::tuple<int, int, double>>> thread_flows(num_threads);
    std::vector<std::vector<double>> thread_times(num_threads);
    std::vector<Eigen::VectorXd> threads_rhs(num_threads, Eigen::VectorXd::Zero(n));


    for (int t = 0; t < this->iteration_count; ++t) {


        auto start = std::chrono::high_resolution_clock::now();

        #pragma omp parallel for schedule(static, chunk)
        for (int u = 0; u < n; u ++) {

#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif
            auto &local_times = thread_times[tid];
            auto &local_flows = thread_flows[tid];
            auto &local_rhs = threads_rhs[tid];

            local_rhs.setZero(n);
            local_rhs[u] =  1.0;
            local_rhs[x_fixed] = -1.0;

            // solve this batch (manual multi; still reuses AMG hierarchy)
            auto t0 = std::chrono::high_resolution_clock::now();
            auto* amgMT = static_cast<AMGSolverMT*>(amg.get());
            const auto& x = (amgMT)->solve_single_thread(local_rhs);
            auto t1 = std::chrono::high_resolution_clock::now();

            const double avg_ms =
                std::chrono::duration<double, std::milli>(t1 - t0).count() / local_rhs.size();
            // store per-RHS times thread-locally (no locks)
            local_times.insert(local_times.end(), local_rhs.size(), avg_ms);

            // compute flows for every result; keep it thread-local
            // (SoA edge list is faster; this uses your existing edges[e].first/second)
            for (int e = 0; e < m; ++e) {
                const int a = edges[e].first;
                const int b = edges[e].second;
                const double fval = edge_weights[e] * (x[a] - x[b]);
                if (std::abs(fval) > 1e-32)
                    local_flows.emplace_back(e, u, fval);
            }

            local_rhs.setZero(n);
        } // end omp parallel for over batches



        // --- merge phase ---
        // Merge all thread buffers (lock-free, outside parallel region)
        for (int tid = 0; tid < num_threads; ++tid) {
            for (double tval : thread_times[tid])
                pure_oracle_running_times.push_back(tval);

            for (auto &[e, u, fval] : thread_flows[tid])
                addFlow(e, u, fval);
        }

        for (auto& vec : thread_times) {
            pure_oracle_running_times.insert(
                pure_oracle_running_times.end(), vec.begin(), vec.end()
            );
        }

        // (Optional) debug: print flows per (edge, source)
        if (debug) {
            printFlow();
        }

        getApproxLoad(load);
        updateEdgeDistances(load);

        auto end = std::chrono::high_resolution_clock::now();
        oracle_running_times.push_back(
            std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
        );

    }


    scaleFlowDown();
    // enforceConservation();
    if (debug) printFlow();
}


void ElectricalFlowParallelOnTheFly::scaleFlowDown() {

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

void ElectricalFlowParallelOnTheFly::getApproxLoad(std::vector<double>& load) {
    const int ell = X.cols();
    load.assign(m, 0.0);

    // --- Buffer for per-edge differences ---
    edge_diffs.resize(static_cast<size_t>(m) * ell);

    // 1. Try dynamic cast to parallel AMG solver
    AMGSolverMT* amg_mt = dynamic_cast<AMGSolverMT*>(amg.get());
    Eigen::MatrixXd potentials;

    if (amg_mt == nullptr) {
        throw std::runtime_error("[ElectricalFlowParallel::getApproxLoad] AMG solver is not parallel-capable (AMGSolverMT).");
    }

    // ---- (A) Parallel batch solve using AMG pool ----
    if (debug)
        std::cout << "[ElectricalFlowParallel] Using parallel AMG batch solver (" << ell << " RHS)\n";

    Eigen::MatrixXd RHS = X; // n × ell (already sparseView->dense fine)
    auto start = std::chrono::high_resolution_clock::now();
    potentials = amg_mt->solve_many(RHS); // parallel multi-RHS solve
    auto end = std::chrono::high_resolution_clock::now();
    pure_oracle_running_times.push_back(
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
    );

    // 2. Compute edge differences in parallel
    const double* __restrict psol = potentials.data();
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int e = 0; e < m; ++e) {
        const int u = edges[e].first;
        const int v = edges[e].second;
        for (int i = 0; i < ell; ++i) {
            edge_diffs[static_cast<size_t>(e) * ell + i] =
                std::fabs(psol[u + i * n] - psol[v + i * n]);
        }
    }

    // 3. Median aggregation and load computation
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int e = 0; e < m; ++e) {
        double* __restrict arr = &edge_diffs[static_cast<size_t>(e) * ell];
        std::nth_element(arr, arr + (ell >> 1), arr + ell);
        const double med = arr[ell >> 1];
        load[e] = edge_weights[e] * med;
    }

    if (debug) {
        double avg_load = std::accumulate(load.begin(), load.end(), 0.0) / m;
        std::cout << "[ElectricalFlowParallel] Avg edge load = " << avg_load << "\n";
    }
}





void ElectricalFlowParallelOnTheFly::storeFlow(){
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



double ElectricalFlowParallelOnTheFly::getFlowForCommodity(int edge_id, int source, int target) {



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

Eigen::SparseMatrix<double> ElectricalFlowParallelOnTheFly::getSketchMatrix(int _m, int _n, double eps) {
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