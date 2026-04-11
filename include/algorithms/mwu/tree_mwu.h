//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_TREE_MWU_H
#define OBLIVIOUSROUTING_TREE_MWU_H

#include "../solver.h"
#include "mwu_framework.h"
#include "oracle/tree/tree_oracle.h"
#include "oracle/tree/tree_transform.h"

// Enum for cycle removal strategy
enum class CycleRemovalStrategy {
    NONE,           // No cycle removal
    TARJAN_SCC,     // Remove cycles using Tarjan SCC algorithm
    NAIVE          // Remove cycles linearly
};

const static std::map<CycleRemovalStrategy, std::string> map_cycle_strategy{
    {CycleRemovalStrategy::NONE, "none"},
    {CycleRemovalStrategy::TARJAN_SCC, "Tarjan_SCC"},
    {CycleRemovalStrategy::NAIVE, "naive"},
};

/*
 *
 * This is the implementation of the tree-based oblivious routing algorithm using the Multiplicative Weights Update (MWU) framework presented by Räcke in 2008.
 * The algorithm iteratively computes a tree (FlatHST) embedding of the graph using the provided oracle, computes the load on the edges based on the tree, and updates the edge distances multiplicatively
 * until the total weight (lambda_sum) of the trees added to the routing table reaches 1.
 */
template<typename HSTDatastructures>
class TreeMWU : public LinearObliviousSolverBase, public MWUFramework {

    std::unique_ptr<TreeOracle<HSTDatastructures>> oracle;
    TreeTransform transform;
    CycleRemovalStrategy cycle_strategy;
    double lambda_sum;
    std::vector<double> rload_current, rload_total;
    std::vector<double> current_distances;
    std::vector<double> mendel_scaling_times;
    std::vector<int> scales;
    std::vector<int> tree_heights;
    double cycleCancellation_time;

public:
    TreeMWU(IGraph& g, int root, std::unique_ptr<TreeOracle<HSTDatastructures>> _oracle,
            CycleRemovalStrategy strategy = CycleRemovalStrategy::NAIVE)
        : LinearObliviousSolverBase(g, root)
        , oracle(std::move(_oracle))
        , transform(graph)
        , cycle_strategy(strategy) {
        current_distances.assign(graph.getNumDirectedEdges(), 1.0);
        rload_current.assign(graph.getNumDirectedEdges(), 0.0);
        rload_total.assign(graph.getNumDirectedEdges(), 0.0);
        lambda_sum = 0.0;
        cycleCancellation_time = 0.0;
    }

    void computeBasisFlows(LinearRoutingTable& table) override {
        table.init(graph);
        run(table);
    }

    void printAdditionalStats() override {
        double average_tree_height = 0.0;
        for (int h : tree_heights) average_tree_height += h;

        std::cout << "Partitioning stats" << std::endl;
        std::cout << "Number of PQ pushes: " << oracle->number_pq_pushes << std::endl;
        std::cout << "Number of PQ pops: " << oracle->number_pq_pushes << std::endl;
        std::cout << "Number of successful relaxations: " << oracle->number_successful_relaxations << std::endl;
        std::cout << "Number of touched nodes: " << oracle->number_touched_nodes << std::endl;

        std::cout << "Average tree height: " << average_tree_height/tree_heights.size() << "\n";
        auto it = ::map_cycle_strategy.find(cycle_strategy);
        if (it != ::map_cycle_strategy.end()) {
            std::cout << "Using cycle cancellation strategy: " << it->second << std::endl;
        }
        std::cout << "Cycle cancellation time: " << this->getCycleCancellationTime() << " micro seconds\n";
        if (oracle->applyMendelScaling) {
            double total = 0.0;
            for (auto t : mendel_scaling_times) total += t;
            std::cout << "Total time spent on Mendel scaling: " << total << "  micro seconds\n";
            double avg = mendel_scaling_times.empty() ? 0.0 : total / mendel_scaling_times.size();
            std::cout << "Average time spent on Mendel scaling per iteration: " << avg << "  micro seconds\n";
        }
    }

    void run(LinearRoutingTable& table) {

        lambda_sum = 0.0;
        while (lambda_sum < 1.0) {
            lambda_sum += treeOracle(table);
            iteration_count++;
        }

        auto t0 = timeNow();
/*
        // Apply cycle removal based on selected strategy
        switch (cycle_strategy) {
            case CycleRemovalStrategy::NONE:
                // No cycle removal
                break;
            case CycleRemovalStrategy::TARJAN_SCC:
                transform.removeCyclesUsingTarjanSCCAlgorithm(table);
                break;
            case CycleRemovalStrategy::NAIVE:
                transform.removeCyclesLinear(table);
                break;
        }
*/
        cycleCancellation_time += duration(timeNow() - t0);
        std::cout << "Path -based table:" << std::endl;

    }

    double treeOracle(LinearRoutingTable& table) {
        auto t0 = timeNow();

        HSTDatastructures t = oracle->getTree(current_distances);
        double oracle_time = duration(timeNow() - t0);
        this->oracle_running_times.push_back(oracle_time);

        int height = 0;
        if constexpr (std::is_same_v<HSTDatastructures, std::shared_ptr<HSTNode>>) {
            height = calculatePointerHSTHeight(t);
        } else {
            height = calculateFlatHSTHeight(t);
        }

        tree_heights.push_back(height);
        scales.push_back(oracle->scales.empty() ? 0.0 : oracle->scales.size());

        t0 = timeNow();
        computeRLoads(t);
        double l = getMaxRload();
        this->load_computation_time += duration(timeNow() - t0);

        t0 = timeNow();
        double lambda = std::min(1.0 / l, 1.0 - lambda_sum);
        computeNewDistances(lambda);
        double weight_update_time = duration(timeNow() - t0);
        this->mwu_weight_update_time += weight_update_time;

        // solve_time = oracle_time + compute time (computeRLoads + computeNewDistances)
        solve_time += oracle_time;

        t0 = timeNow();

        TreeIteration<HSTDatastructures> iter(std::move(t), current_distances, lambda);
        transform.transform(iter, table);
        this->transformation_time += duration(timeNow() - t0);

        if (oracle->applyMendelScaling)
            mendel_scaling_times.push_back(oracle->getMendelScalingTime());

        return lambda;
    }

    void computeRLoads(const FlatHST& hst) {
        std::fill(rload_current.begin(), rload_current.end(), 0.0);

        std::queue<int> q;
        q.push(hst.root());
        while (!q.empty()) {
            int node_idx = q.front(); q.pop();
            for (int child_idx : hst.children(node_idx)) {
                q.push(child_idx);

                auto node_members  = hst.memberRange(node_idx);
                auto child_members = hst.memberRange(child_idx);
                if (child_members.size() == node_members.size() || child_members.empty()) continue;

                std::vector<char> S(graph.getNumNodes(), 0);
                for (int v : child_members) S[v] = 1;

                double cut = 0.0;
                for (int u : child_members)
                    for (int v : graph.neighbors(u))
                        if (!S[v]) cut += graph.getEdgeCapacity(u, v);
                if (cut <= 1e-12) cut = 1e-12;

                int repParent = node_members.empty() ? child_members[0] : node_members[0];
                int repChild  = child_members[0];

                auto path = graph.getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int u = path[i], v = path[i+1];
                    int e = graph.getEdgeId(u, v), anti_e = graph.getEdgeId(v, u);
                    double cap = std::max(graph.getEdgeCapacity(u, v), 1e-12);
                    double rload = rload_current[e] + cut / cap;
                    rload_current[e] = rload_current[anti_e] = rload;
                }
            }
        }
    }

    void computeRLoads(const std::shared_ptr<HSTNode>& t) {
        assert(t != nullptr);
        std::fill(rload_current.begin(), rload_current.end(), 0.0);

        std::queue<const HSTNode*> q;
        q.push(t.get());
        while (!q.empty()) {
            const HSTNode* node = q.front(); q.pop();
            for (const auto& child_ptr : node->getChildren()) {
                const HSTNode* child = child_ptr.get();
                q.push(child);

                if (node->getMembers().size() == child->getMembers().size() || child->getMembers().empty()) continue;
                const std::vector<int>& clusterVertices = child->getMembers();

                std::vector<char> S(graph.getNumNodes(), 0);
                for (int v : clusterVertices) S[v] = 1;

                double cut = 0.0;
                for (int u : clusterVertices)
                    for (auto v : graph.neighbors(u))
                        if (!S[v]) cut += graph.getEdgeCapacity(u, v);
                if (cut <= 1e-12) cut = 1e-12;

                int repParent = node->getMembers().empty() ? clusterVertices[0] : node->getMembers()[0];
                int repChild  = clusterVertices[0];

                auto path = graph.getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int u = path[i], v = path[i+1];
                    int e = graph.getEdgeId(u, v), anti_e = graph.getEdgeId(v, u);
                    double cap = std::max(graph.getEdgeCapacity(u, v), 1e-12);
                    double rload = rload_current[e] + cut / cap;
                    rload_current[e] = rload_current[anti_e] = rload;
                }
            }
        }
    }

    double getMaxRload() const {
        assert(rload_current.size() > 0);

        double max_load = 0.0;
        for (const auto& load : rload_current) {
            if (load > max_load) {
                max_load = load;
            }
        }
        return max_load;
    }

    void addCurrentLoad(double lambda) {
        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
            rload_total[e] += rload_current[e] * lambda;
        }
    }
    void computeNewDistances(double lambda) {
        addCurrentLoad(lambda);

        double max_r = 0.0;
        for (auto& r : rload_total) max_r = std::max(max_r, r);

        double sumExp = 0.0;
        for (auto& r : rload_total) sumExp += std::exp(r - max_r);
        if (sumExp <= 0.0 || !std::isfinite(sumExp)) sumExp = 1.0;

        double min_d = std::numeric_limits<double>::infinity();
        std::vector<double> newDist;
        newDist.resize(rload_total.size());

        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
            double r = rload_total[e];

            double cap = graph.getEdgeCapacity(e);
            if (cap < EPS) cap = EPS;
            double d = (std::exp(r - max_r) / cap) / sumExp;
            newDist[e] = d;
            if (d < min_d) min_d = d;
        }
        if (min_d < EPS) min_d = EPS;

        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
            double d = newDist[e];
            double norm = d / min_d;
            if (norm < 1.0) norm = 1.0;
            current_distances[e] = norm;
        }

        updateDistances(current_distances);
    }

    void updateDistances(const std::vector<double> &distances) override {
        for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
            graph.updateEdgeDistance(e, distances[e]);
        }
    }

    std::vector<int> getScales() const {
        return scales;
    }

    double getCycleCancellationTime() const {
        return cycleCancellation_time;
    }
};

#endif //OBLIVIOUSROUTING_TREE_MWU_H

