//
// Created by Mert Biyikli on 04.02.26.
//

#include "electrical_flow_par_batches.h"
#include <omp.h>
#include "../../utils/my_math.h"
#include <random>
#include <thread>

void ElectricalFlowParallelBatches::init(bool debug)
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


    initEdgeDistances();
    buildWeightDiag();
    buildIncidence();

    // init AMG
    int thread_nums =0;
    #ifdef _OPENMP
    thread_nums = omp_get_max_threads();
    #endif
    p_threads = std::max(1, thread_nums);

    amg_pool.clear();
    amg_pool.reserve(p_threads);

    for (int t = 0; t < p_threads; ++t) {
        auto s = std::make_unique<AMGSolver>();
        s->init(graph, edge_weights, n, edges, debug);
        amg_pool.emplace_back(std::move(s));
    }

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

void ElectricalFlowParallelBatches::initEdgeDistances() {
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
        double cap = graph.getEdgeCapacity(e);      // undirected capacity accessor
        edge_capacities[e]  = cap;
        edge_probabilities[e] = edge_distances[e] / cap_X;
        edge_weights[e] = std::pow(cap, 2) / (edge_probabilities[e] + inv_m);
    }


}


void ElectricalFlowParallelBatches::extract_edge_list() {
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



void ElectricalFlowParallelBatches::buildIncidence()
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

void ElectricalFlowParallelBatches::buildWeightDiag() {
    if (!W.nonZeros()) {
        W = Eigen::SparseMatrix<double>(m, m);
        std::vector<Eigen::Triplet<double>> w_triplets; w_triplets.reserve(m);
        for (int e = 0; e < m; ++e)
            w_triplets.emplace_back(e, e, edge_weights[e]);
        W.setFromTriplets(w_triplets.begin(), w_triplets.end());
    }
}

void ElectricalFlowParallelBatches::refreshWeightMatrix() {
    if (!W.nonZeros()) {
        buildWeightDiag();
    } else {
        for (int e = 0; e < m; ++e)
            W.coeffRef(e, e) = edge_weights[e];
    }
}

void ElectricalFlowParallelBatches::updateEdgeDistances(const std::vector<double>& load) {
    // --- 1) stable MWU in log-space with damping & clipping ---
    const double eta   = 1.0 / (4.0 * roh);   // smaller than your 1/(2*roh)
    const double gclip = 10.0;                // gradient clip (tune 5–20)

    // Kahan sum for cap_X
    double capX = 0.0, c = 0.0;

    #ifdef _OPENMP
    #pragma omp parallel for reduction(+:cap_X) schedule(dynamic)
    #endif
    for (int e = 0; e < m; ++e) {
        if (load[e] <= 0.0 || std::isnan(load[e])) continue;
        edge_distances[e] *= (1 + (1/(2*roh)) * load[e]);
        cap_X += edge_distances[e];
    }

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


    for (auto &solver : amg_pool) {
        solver->updateAllEdges(edge_weights, edges);
        solver->updateSolver();
    }


    refreshWeightMatrix();
}




void ElectricalFlowParallelBatches::run(LinearRoutingTable &table) {
    // Sparse RHS with only two nonzeros per demand (+1 at u, -1 at x_fixed)
    std::vector<double> load(m, 0.0);


    std::vector<std::vector<double>> div_local(p_threads, std::vector<double>(n,0.0));

    // each thread has its own flow table
    std::vector<LinearRoutingTable> thread_flows(p_threads);
    for (auto& tf : thread_flows) {
        tf.init(graph);
    }
    std::vector<double> thread_times(p_threads);

    auto t0 = timeNow();

    for (int t = 0; t < this->iteration_count; ++t) {
        auto start_time = std::chrono::high_resolution_clock::now();
        for (int tid = 0; tid < p_threads; ++tid) {
            for ( auto& v : thread_flows[tid].src_flows ) {
                v.clear();
            }
            for ( auto& v : thread_flows[tid].src_ids ) {
                v.clear();
            }
            std::fill(div_local[tid].begin(), div_local[tid].end(), 0.0);
        }
        std::fill(div_accum.begin(), div_accum.end(), 0.0);

        // --- main loop over sources (u -> x_fixed) ---
        // --- main loop over sources (u -> x_fixed) ---
        // NO OpenMP here. We do outer std::thread parallelism and force AMGCL to 1 thread inside each worker.

        const int p = p_threads; // you set this in init()
        std::vector<std::thread> workers;
        workers.reserve(p);

        auto worker_fn = [&](int tid) {
            Eigen::setNbThreads(1);     // avoid Eigen oversubscription too

            auto &local_flows = thread_flows[tid];
            auto &div         = div_local[tid];


            const int start = (tid * n) / p;
            const int end   = ((tid + 1) * n) / p;

            Eigen::VectorXd rhs(n);
            Eigen::VectorXd pot(n);

            for (int u = start; u < end; ++u) {
                if (u == x_fixed) continue;

                rhs.setZero();
                rhs(u) = 1.0;
                rhs(x_fixed) = 0.0; // you ground x_fixed by your Dirichlet procedure / matrix modification

                // solve L * pot = rhs
                pot = (*amg_pool[tid]).solve(rhs, epsilon_L);   // adapt signature if yours differs

                // potentials -> flows
                for (int e = 0; e < m; ++e) {
                    const int a = edges[e].first;
                    const int b = edges[e].second;

                    const double fval = edge_weights[e] * (pot[a] - pot[b]);
                    if (std::abs(fval) < EPS) continue;

                    if (u == b)      div[u] += fval;
                    else if (u == a) div[u] -= fval;

                    const int id   = graph.getEdgeId(a, b);
                    const int anti = graph.getAntiEdge(id);

                    if (fval >= 0) local_flows.addFlow(id,   u, std::abs(fval));
                    else           local_flows.addFlow(anti, u, std::abs(fval));
                }
            }
        };

        for (int tid = 0; tid < p; ++tid) workers.emplace_back(worker_fn, tid);
        for (auto &th : workers) th.join();

        // --- merge phase ---
        // Merge all thread buffers (lock-free, outside parallel region)
        double oracle_time = duration(timeNow() - start_time);
        oracle_running_times.push_back(oracle_time);


        for (int tid = 0; tid < p; ++tid) {
            auto& local_table = thread_flows[tid];
            for (int e = 0; e < local_table.src_ids.size(); ++e) {
                for (int j = 0; j < local_table.src_ids[e].size(); ++j) {
                    int source = local_table.src_ids[e][j];
                    double flow = local_table.src_flows[e][j];
                    if (std::abs(flow)>EPS) {
                        table.addFlow(e, source, flow);
                    }
                }
            }
        }

            for (auto div_thread  :div_local) {
                for (int i = 0; i < n; ++i) {
                    div_accum[i] += div_thread[i];
                }
            }



        getApproxLoad(load);
        updateEdgeDistances(load); // keeps 'weights' vector in sync

        oracle_running_times.push_back(duration(timeNow() - start_time));

    }

    this->solve_time = duration(timeNow() - t0);
}

void ElectricalFlowParallelBatches::scaleFlowDown(LinearRoutingTable& table) {
    // scale the flow from the adjacency list flow
    if (iteration_count > 0) {
        const double inv_iters = 1.0 / static_cast<double>(iteration_count);
        for (int e = 0; e < graph.getNumEdges(); ++e) // dont use m here. m is undirected edges only
            for (double &val : table.src_flows[e]) val *= inv_iters;
    }
}

void ElectricalFlowParallelBatches::getApproxLoad(std::vector<double>& load) {
    const int ell = X.cols();

    // 0) epsilon_L used by solver
    epsilon_L = epsilon / (8.0 * edges.size() * std::pow(n, 4) * K);

    // 1) Build RHS (n x ell)
    Eigen::MatrixXd RHS = Eigen::MatrixXd(X);   // X is sparse; you were already densifying implicitly
    RHS.row(x_fixed).setZero();

    // 2) Solve potentials (n x ell) in parallel over RHS columns
    Eigen::MatrixXd potentials(n, ell);
    potentials.setZero();

    const int p = p_threads;
    std::vector<std::thread> workers;
    workers.reserve(p);

    auto solve_worker = [&](int tid) {
    #ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(1);
    #endif
        Eigen::setNbThreads(1);


        const int c0 = (tid * ell) / p;
        const int c1 = ((tid + 1) * ell) / p;

        Eigen::VectorXd rhs(n), sol(n);

        for (int j = c0; j < c1; ++j) {
            rhs = RHS.col(j);
            sol = (*amg_pool[tid]).solve(rhs, epsilon_L);     // adapt to your exact signature
            potentials.col(j) = sol;
        }
    };

    for (int tid = 0; tid < p; ++tid) workers.emplace_back(solve_worker, tid);
    for (auto& th : workers) th.join();

    // 3) edge_diffs[e,i] = |pot[u,i] - pot[v,i]|
    edge_diffs.resize(static_cast<size_t>(m) * ell);

    auto diffs_worker = [&](int tid) {
    #ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(1);
    #endif
        Eigen::setNbThreads(1);

        const int e0 = (tid * m) / p;
        const int e1 = ((tid + 1) * m) / p;

        for (int e = e0; e < e1; ++e) {
            const int u = edges[e].first;
            const int v = edges[e].second;

            double* __restrict out = &edge_diffs[static_cast<size_t>(e) * ell];
            for (int i = 0; i < ell; ++i) {
                out[i] = std::fabs(potentials(u, i) - potentials(v, i));
            }
        }
    };

    workers.clear();
    for (int tid = 0; tid < p; ++tid) workers.emplace_back(diffs_worker, tid);
    for (auto& th : workers) th.join();

    // 4) median per edge + load
    auto median_worker = [&](int tid) {
    #ifdef _OPENMP
        omp_set_dynamic(0);
        omp_set_num_threads(1);
    #endif
        Eigen::setNbThreads(1);

        const int e0 = (tid * m) / p;
        const int e1 = ((tid + 1) * m) / p;

        for (int e = e0; e < e1; ++e) {
            double* __restrict arr = &edge_diffs[static_cast<size_t>(e) * ell];
            std::nth_element(arr, arr + (ell >> 1), arr + ell);
            const double med = arr[ell >> 1];
            load[e] = edge_weights[e] * med;
        }
    };

    workers.clear();
    for (int tid = 0; tid < p; ++tid) workers.emplace_back(median_worker, tid);
    for (auto& th : workers) th.join();
}




// Compute the median absolute difference between two node potentials
// across all ℓ sampled solutions.
//   diffs_u = potentials[u][*]  (length ℓ)
//   diffs_v = potentials[v][*]  (length ℓ)
/*
double ElectricalFlowParallelBatches::recoverNorm(const std::vector<double>& diffs_u,
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
*/

Eigen::SparseMatrix<double> ElectricalFlowParallelBatches::getSketchMatrix(int _m, int _n, double eps) {
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