#include "../../mwu/electrical_mwu.h"

template <typename ExecutionPolicy>
class ParElectricalFlowMWU : public ElectricalMWU {
    ExecutionPolicy executionStrategy;
public:
    ParElectricalFlowMWU(IGraph& g, int root, bool use_sketching, ExecutionPolicy _executionStrategy)
    : ElectricalMWU(g, root, use_sketching, true), executionStrategy(std::move(_executionStrategy)) {
    }

    std::mutex mutex;

    void run(LinearRoutingTable &table) override {
        auto t0 = timeNow();
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
        Eigen::VectorXd potentials(n);
        std::vector<double> load(m, 0.0);
        solve_time += duration(timeNow() - t0);


        for (int t = 0; t < this->iteration_count; ++t) {

            auto oracle_iteration = 0.0;
            executionStrategy.forEachSource(graph.getNumNodes(), [&](int source) {
                auto local_result = sourceRouting(source, rhs, potentials);

                std::lock_guard<std::mutex> lock(mutex);

            });
            postRun(load, oracle_iteration);
        }
    }

    LinearRoutingTable sourceRouting(int source, Eigen::VectorXd& rhs, Eigen::VectorXd& potentials) {
        LinearRoutingTable local_table;
        local_table.init(graph);
        rhs = buildDemandVector(source);

        laplaceSolve(rhs, potentials);
        // addFlowToTable measures its own time and adds to transformation_time
        addFlowToTable(source, potentials, local_table);

        return local_table;
    }

    Eigen::VectorXd buildDemandVector(int source) {
        Eigen::VectorXd demandVector(n);
        demandVector.setZero();
        demandVector[source] = 1.0;
        demandVector[x_fixed] = -1.0;
        return demandVector;
    }

    void laplaceSolve(const Eigen::VectorXd& rhs, Eigen::VectorXd& potentials) {
        auto t0 = timeNow();
        potentials = amg->solve(rhs);
        solve_time += duration(timeNow() - t0);
    }

    void postRun(std::vector<double>& load, double oracle_iteration) {
        auto t0 = timeNow();
        if (use_sketching) {
            getApproxLoad(load);
        }else {
            getExactLoad(load);
        }
        this->load_computation_time += duration(timeNow()-t0);

        t0 = timeNow();
        updateDistances(load);
        double weight_update_time = duration(timeNow() - t0);
        mwu_weight_update_time += weight_update_time;

        oracle_running_times.push_back(oracle_iteration);
    }

};



struct SequentialExecution {
    template <typename Task>
    void forEachSource(int numSources, Task&& task) const {
        for (int s = 0; s < numSources; ++s) {
            task(s);
        }
    }
};


struct OpenMPExecution {
    template <typename Task>
    void forEachSource(int numSources, Task&& task) const {
        #pragma omp parallel for schedule(dynamic)
        for (int s = 0; s < numSources; ++s) {
            task(s);
        }
    }
};