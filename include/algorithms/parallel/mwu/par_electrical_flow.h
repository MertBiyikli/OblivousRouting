#pragma once

#include "../../mwu/electrical_mwu.h"
#include "../../../utils/my_math.h"

#include <vector>
#include <memory>
#include <numeric>
#include <algorithm>
#include <cmath>

template <typename ExecutionPolicy>
class ParElectricalFlowMWU final : public ElectricalMWU {
private:
    ExecutionPolicy execution;

    std::vector<std::unique_ptr<LaplacianSolver>> workerSolvers;
    std::vector<LinearRoutingTable> workerTables;

public:
    ParElectricalFlowMWU(
        IGraph& g,
        int root,
        bool useSketching,
        ExecutionPolicy executionPolicy
    )
        : ElectricalMWU(g, root, useSketching, false),
          execution(std::move(executionPolicy)) {}

    void run(LinearRoutingTable& table) override {
        configureOuterOpenMP(execution.numWorkers());
        LaplacianSolver::printOpenMPDiagnostics("ParElectricalFlowMWU::run before initWorkerTables");

        initWorkerTables();

        std::vector<double> load(m, 0.0);

        for (int iter = 0; iter < iteration_count; ++iter) {
            const double oracleTime = routeAllSourcesParallel();

            computeLoadParallel(load);

            updateDistancesParallel(load);

            updateAllWorkerSolvers();

            oracle_running_times.push_back(oracleTime);
        }

        mergeWorkerTables(table);
    }

    void initAMGSolver(boost::property_tree::ptree params) override {
        const int workers = execution.numWorkers();

        workerSolvers.clear();
        workerSolvers.reserve(workers);

        for (int w = 0; w < workers; ++w) {
            auto solver = std::make_unique<LaplacianSolver>();
            solver->setSolverParams(params);
            solver->init(graph, edge_weights, n, edges, debug);
            workerSolvers.push_back(std::move(solver));
        }
    }

private:
    void initWorkerTables() {
        const int workers = execution.numWorkers();

        workerTables.assign(workers, LinearRoutingTable{});

        for (auto& localTable : workerTables) {
            localTable.init(graph);
        }
    }

    double routeAllSourcesParallel() {
        const int workers = execution.numWorkers();

        std::vector<double> localSolveTimes(workers, 0.0);
        std::vector<double> localTransformTimes(workers, 0.0);

        const auto wallStart = timeNow();

        execution.parallelFor(0, n, [&](int source, int workerId) {
            assert(workerId >= 0);
            assert(workerId < workers);
            assert(workerId < static_cast<int>(workerSolvers.size()));
            assert(workerId < static_cast<int>(workerTables.size()));

            if (source == x_fixed) {
                return;
            }

            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
            rhs[source] = 1.0;
            rhs[x_fixed] = -1.0;

            auto t0 = timeNow();
            Eigen::VectorXd potentials = workerSolvers[workerId]->solve(rhs);
            const double solveDuration = duration(timeNow() - t0);

            localSolveTimes[workerId] += solveDuration;

            t0 = timeNow();
            addFlowToTableNoTiming(source, potentials, workerTables[workerId]);
            localTransformTimes[workerId] += duration(timeNow() - t0);
        });

        const double wallOracleTime = duration(timeNow() - wallStart);

        double cpuSolveTime = 0.0;
        for (double value : localSolveTimes) {
            cpuSolveTime += value;
        }

        double cpuTransformTime = 0.0;
        for (double value : localTransformTimes) {
            cpuTransformTime += value;
        }

        // For benchmark comparison, use wall time.
        solve_time += wallOracleTime;

        // Transformation is still accumulated as worker CPU time for now.
        transformation_time += cpuTransformTime;

        return wallOracleTime;
    }

    void addFlowToTableNoTiming(
        int source,
        const Eigen::VectorXd& potential,
        LinearRoutingTable& table
    ) {
        for (int e = 0; e < m; ++e) {
            const auto& [head, tail] = edges[e];

            const double fval =
                edge_weights[e] *
                (potential[head] - potential[tail]);

            if (std::abs(fval) > SOFT_EPS) {
                const int originalEdgeId = graph.getEdgeId(head, tail);

                const int directedEdge =
                    fval < 0.0
                    ? graph.getAntiEdge(originalEdgeId)
                    : originalEdgeId;

                table.addFlow(directedEdge, source, std::abs(fval));
            }
        }
    }

    void computeLoadParallel(std::vector<double>& load) {
        auto t0 = timeNow();

        if (use_sketching) {
            getApproxLoadParallel(load);
        } else {
            getExactLoadParallel(load);
        }

        load_computation_time += duration(timeNow() - t0);
    }

    void getExactLoadParallel(std::vector<double>& load) {

        const int workers = execution.numWorkers();

        std::vector<std::vector<double>> localLoads(
            workers,
            std::vector<double>(m, 0.0)
        );

        std::vector<double> localSolveTimes(workers, 0.0);

        execution.parallelFor(0, m, [&](int f, int workerId) {
            assert(workerId >= 0);
    assert(workerId < workers);
    assert(workerId < static_cast<int>(workerSolvers.size()));
    assert(workerId < static_cast<int>(localLoads.size()));
    assert(workerId < static_cast<int>(localSolveTimes.size()));

            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);

            const auto& [u_f, v_f] = edges[f];
            rhs[u_f] = 1.0;
            rhs[v_f] = -1.0;

            auto t0 = timeNow();
            Eigen::VectorXd sol = workerSolvers[workerId]->solve(rhs);
            localSolveTimes[workerId] += duration(timeNow() - t0);

            auto& localLoad = localLoads[workerId];

            for (int e = 0; e < m; ++e) {
                const auto& [u_e, v_e] = edges[e];
                const double potDiff = sol[v_e] - sol[u_e];

                localLoad[e] += edge_weights[e] * std::abs(potDiff);
            }
        });

        std::fill(load.begin(), load.end(), 0.0);

        for (int worker = 0; worker < workers; ++worker) {
            for (int e = 0; e < m; ++e) {
                load[e] += localLoads[worker][e];
            }

            solve_time += localSolveTimes[worker];
        }
    }

    void getApproxLoadParallel(std::vector<double>& load) {
        const int ell = X.cols();
        const int workers = execution.numWorkers();

        edge_diffs.assign(static_cast<size_t>(m) * ell, 0.0);

        std::vector<double> localSolveTimes(workers, 0.0);

        execution.parallelFor(0, ell, [&](int i, int workerId) {
            assert(workerId >= 0);
            assert(workerId < workers);
            assert(workerId < static_cast<int>(workerSolvers.size()));
            assert(workerId < static_cast<int>(localSolveTimes.size()));

            Eigen::VectorXd rhs = X.col(i);
            auto t0 = timeNow();
            Eigen::VectorXd sol = workerSolvers[workerId]->solve(rhs);
            localSolveTimes[workerId] += duration(timeNow() - t0);

            for (int e = 0; e < m; ++e) {
                const auto& [u, v] = edges[e];
                const double diff = sol[v] - sol[u];

                edge_diffs[static_cast<size_t>(e) * ell + i] =
                    std::abs(diff);
            }
        });

        for (double value : localSolveTimes) {
            solve_time += value;
        }

        execution.parallelFor(0, m, [&](int e, int workerId) {
            double* arr = &edge_diffs[static_cast<size_t>(e) * ell];

            std::nth_element(
                arr,
                arr + (ell >> 1),
                arr + ell
            );

            const double med = arr[ell >> 1];
            load[e] = edge_weights[e] * med;
        });
    }

    void updateDistancesParallel(const std::vector<double>& load) {
        auto t0 = timeNow();

        const int workers = execution.numWorkers();
        std::vector<double> localCap(workers, 0.0);

        execution.parallelFor(0, m, [&](int e, int workerId) {
            assert(workerId >= 0);
            assert(workerId < static_cast<int>(workerSolvers.size()));
            const double le = load[e];

            if (le > 0.0 && std::isfinite(le)) {
                edge_distances[e] *=
                    1.0 + (1.0 / (2.0 * roh)) * le;
            }

            localCap[workerId] += edge_distances[e];
        });

        cap_X = 0.0;
        for (double value : localCap) {
            cap_X += value;
        }

        execution.parallelFor(0, m, [&](int e, int workerId) {
            const double x = edge_distances[e];

            if (x <= 0.0 || std::isnan(x)) {
                return;
            }

            const double p = x / cap_X;
            edge_probabilities[e] = p;

            const double cap = edge_capacities[e];
            edge_weights[e] = std::pow(cap, 2) / (p + inv_m);
        });

        mwu_weight_update_time += duration(timeNow() - t0);
    }

    void updateAllWorkerSolvers() {

        const int workers = execution.numWorkers();
        execution.parallelFor(
            0,
            static_cast<int>(workerSolvers.size()),
            [&](int solverIndex, int workerId) {
                assert(workerId >= 0);
                assert(workerId < workers);
                assert(solverIndex >= 0);
                assert(solverIndex < static_cast<int>(workerSolvers.size()));

                workerSolvers[solverIndex]->updateAllEdges(edge_weights, edges);
                workerSolvers[solverIndex]->updateSolver();
            }
        );
    }

    void mergeWorkerTables(LinearRoutingTable& mainTable) {
        for (const auto& localTable : workerTables) {
            for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
                for (int source = 0; source < n; ++source) {
                    const double value = localTable.getFlow(e, source);

                    if (std::abs(value) > EPS) {
                        mainTable.addFlow(e, source, value);
                    }
                }
            }
        }
    }

    void configureOuterOpenMP(int threads) {
#ifdef OR_ENABLE_OPENMP
        omp_set_dynamic(0);
        omp_set_max_active_levels(1);
        omp_set_num_threads(threads);
#endif
    }
};