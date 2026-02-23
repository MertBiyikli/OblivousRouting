//
// Created by Mert Biyikli on 04.02.26.
//

#include "electrical_flow_par_batches.h"
#include "../../utils/my_math.h"
#include <thread>


static inline unsigned default_threads() {
    unsigned t = std::thread::hardware_concurrency();
    return t ? t : 4;
}

void ParallelElectricalMWU::initAMGSolver(boost::property_tree::ptree _params)
{

    // init AMG
    int thread_nums = default_threads();


    p_threads = std::max(1, thread_nums);
    std::cout << "Initializing with " << p_threads << " threads for parallel batches.\n";

    amg_pool.clear();
    amg_pool.reserve(p_threads);

    for (int t = 0; t < p_threads; ++t) {
        auto s = std::make_unique<AMGSolver>();
        s->init(graph, edge_weights, n, edges, debug);
        s->setSolverParams(_params);
        amg_pool.emplace_back(std::move(s));
    }
}

void ParallelElectricalMWU::updateEdgeDistances(const std::vector<double>& load) {
    for (int e = 0; e < m; ++e) {
        if (load[e] <= 0.0 || std::isnan(load[e])) continue;
        edge_distances[e] *= (1 + (1/(2*roh)) * load[e]);
        cap_X += edge_distances[e];
    }

    for (int e = 0; e < m; ++e) {
        double x = edge_distances[e];
        if (x <= 0.0 || std::isnan(x)) continue;

        double p = x / cap_X;
        edge_probabilities[e] = p;

        double cap = edge_capacities[e];
        double w   = std::pow(cap, 2) / (p + inv_m);
        edge_weights[e] = w;
    }

    for (auto& amg_ : amg_pool) {
        amg_->updateAllEdges(edge_weights, edges);
        amg_->updateSolver();
    }
}



void ParallelElectricalMWU::run(LinearRoutingTable &table) {
    // Sparse RHS with only two nonzeros per demand (+1 at u, -1 at x_fixed)
    std::vector load(m, 0.0);
    std::vector<std::vector<double>> div_local(p_threads, std::vector<double>(n,0.0));

    // each thread has its own flow table and time table
    std::vector<double> thread_times(p_threads);
    std::vector<LinearRoutingTable> thread_flows(p_threads);
    for (auto& tf : thread_flows) {
        tf.init(graph);
    }

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
        }
        std::vector<std::thread> workers;
        const int p = p_threads; // you set this in init()
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

                addFlowToTable(u, pot, local_flows);
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

        getApproxLoad(load);
        updateEdgeDistances(load); // keeps 'weights' vector in sync

        oracle_running_times.push_back(duration(timeNow() - start_time));

    }

    this->solve_time = duration(timeNow() - t0);
}


void ParallelElectricalMWU::getApproxLoad(std::vector<double>& load) {
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

    // recover norm
    auto median_worker = [&](int tid) {

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

