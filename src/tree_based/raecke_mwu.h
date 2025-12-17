//
// Created by Mert Biyikli on 08.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_LINEAR_OR_H
#define OBLIVIOUSROUTING_RAECKE_LINEAR_OR_H


#include "raecke_oracle_iteration.h"
#include "../solver/solver.h"
#include "../solver/mwu_framework.h"
#include "raecke_oracle.h"
#include "raecke_oracle_iteration.h"
#include "raecke_tree.h"
#include "ckr/optimized/efficient_raecke_ckr_transform.h"
#include "raecke_linear_routing.h"

class RaeckeMWU : public MWUFramework {
public:
    std::unique_ptr<RaeckeOracle> oracle;
    double lambdaSum;

    std::vector<double> rload_current;
    std::vector<double> rload_total;

    std::vector<OracleTreeIteration> iteration;

    std::vector<double> current_distances;

    RaeckeMWU(IGraph& g, int root,
              std::unique_ptr<RaeckeOracle> _oracle):MWUFramework(g, root),
    oracle(std::move(_oracle)),
          lambdaSum(0.0) {

        current_distances.assign(graph.getNumEdges(), 1.0);
        rload_current.assign(graph.getNumEdges(), 0.0);
        rload_total.assign(graph.getNumEdges(), 0.0);
    }

    // Each Racke Variant has to implement its own transformSolution method
    virtual void transformSolution(LinearRoutingTable& table) = 0;

    void computeBasisFlows(LinearRoutingTable &table) override {
        run();
        normalizeLambdas();
        transformSolution(table);
        //table.printFlows(graph);
        assert(table.isValid(graph));
        // print out the linear routing table

    }

    void run() {
        lambdaSum = 0.0;
        while (lambdaSum < 1.0) {
            auto start = std::chrono::high_resolution_clock::now();
            lambdaSum += treeOracle();
            oracle_running_times.push_back(
                std::chrono::duration<double, std::milli>(
                    std::chrono::high_resolution_clock::now() - start
                ).count()
            );
            iteration_count++;
        }
    }

    double treeOracle() {
        auto t = oracle->getTree(current_distances);
        assert(t != nullptr);
        computeRLoads(t);
        double l = getMaxRload();
        double lambda = std::min(1.0/l, 1.0 - lambdaSum);

        OracleTreeIteration iter(lambda, t, current_distances);
        iteration.push_back(iter);

        computeNewDistances(lambda);

        return lambda;
    }

    void normalizeLambdas() {
        for (auto& it : iteration) {
            double lambda = it.getLambda();
            lambda /= lambdaSum;
            it.setLambda(lambda);
        }
    }


    void addCurrentLoad(double lambda) {
        for (int e = 0; e < graph.getNumEdges(); ++e) {
            rload_total[e] += rload_current[e] * lambda;
        }
    }

    void computeRLoads(std::shared_ptr<ITreeNode> t) {
        assert(t != nullptr);
        // clear current rloads

        std::fill(rload_current.begin(), rload_current.end(), 0.0);



        std::queue<std::shared_ptr<ITreeNode>> q;
        q.push(t);
        while (!q.empty()) {
            auto node = q.front();
            q.pop();

            // --- 3️⃣ Process each child: represents a cut S_child | V\S_child ---
            for (auto child : node->getChildren()) {
                // Add child to traversal queue
                q.push(child);


                // if node has the same nodes as child, skip
                if (node->getMembers().size() == child->getMembers().size()) {
                    continue;
                }
                const std::vector<int>& clusterVertices = child->getMembers();
                if (clusterVertices.empty()) continue;

                // Build set for fast lookup
                std::vector<char> S(graph.getNumNodes(), 0);
                for (int v : clusterVertices) S[v] = 1;


                // --- 4️⃣ Compute total cut capacity of this child cluster ---
                double cut = 0.0;
                for (int u : clusterVertices) {
                    for (auto&  v : graph.neighbors(u)) {
                        if (!S[v])  // boundary edge
                            cut += graph.getEdgeCapacity(u, v);
                    }
                }
                if (cut <= 1e-12) cut = 1e-12;  // avoid zero-division

                // --- 5️⃣ Choose representative vertices (simple heuristic) ---
                int repParent = node->getMembers().empty() ? clusterVertices[0] : node->getMembers()[0];
                int repChild  = clusterVertices[0];

                // --- 6️⃣ Compute shortest path between parent and child reps ---
                auto path = graph.getShortestPath(repParent, repChild);
                if (path.size() < 2) continue;

                // --- 7️⃣ Update edge r-loads along that path ---
                for (size_t i = 0; i + 1 < path.size(); ++i) {
                    int u = path[i];
                    int v = path[i + 1];

                    int e = graph.getEdgeId(u, v);
                    double rload = rload_current[e];
                    double cap = graph.getEdgeCapacity(u, v);
                    if (cap <= 1e-12) cap = 1e-12;

                    rload += cut / cap;
                    rload_current[e] = rload;

                    int anti_e = graph.getEdgeId(v, u);
                    rload_current[anti_e] = rload;

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

    void computeNewDistances(double lambda) {
        addCurrentLoad(lambda);

        double max_r = 0.0;
        for (auto& r : rload_total) max_r = std::max(max_r, r);

        double sumExp = 0.0;
        for (auto& r : rload_total) sumExp += std::exp(r - max_r);
        if (sumExp <= 0.0 || !std::isfinite(sumExp)) sumExp = 1.0;

        // 3) distances
        double min_d = std::numeric_limits<double>::infinity();
        std::vector<double> newDist;
        newDist.reserve(rload_total.size());

        for (int e = 0; e < graph.getNumEdges(); ++e) {
            double r = rload_total[e];

            double cap = graph.getEdgeCapacity(e);
            if (cap < EPS) cap = EPS;
            double d = (std::exp(r - max_r) / cap) / sumExp;
            newDist[e] = d;
            if (d < min_d) min_d = d;
        }
        if (min_d < EPS) min_d = EPS;

        for (int e = 0; e < graph.getNumEdges(); ++e) {
            double d = newDist[e];
            double norm = d / min_d;
            if (norm < 1.0) norm = 1.0;
            graph.updateEdgeDistance(e, norm);
            current_distances[e] = norm;
        }
    }

    void addTree(const OracleTreeIteration& iter, LinearRoutingTable& table) {
        distributeDemands(iter.getTree(), iter, table);
        removeCycles(table);
    }

    void distributeDemands(const std::shared_ptr<ITreeNode> node, const OracleTreeIteration& iter, LinearRoutingTable& table) {
        if (!node) return;

        const double& lambda = iter.getLambda();
        const auto& distance = iter.getDistance();
        for (const auto& child : node->getChildren()) {
            std::set<int> A = collectSubtreeVertices(child);
            std::set<int> B;
            for (int v : graph.getVertices()) {
                if (!A.count(v)) B.insert(v);
            }

            for (int src : A) {
                for (int dst : B) {
                    if (src == dst) continue;
                    auto path = graph.getShortestPath(src, dst, iter.getDistance());

                    pushFlowAlongPath(path, src, dst, lambda, table);
                }
            }

            distributeDemands(child, iter, table);
        }
    }

    void pushFlowAlongPath(const std::vector<int>& p,
                           int s,
                           int t,
                           double lambda,
                           LinearRoutingTable& table) {
        for (size_t i = 0; i + 1 < p.size(); ++i) {
            int u = p[i];
            int v = p[i + 1];

            // first, identify the orientation of the flow
            // we have to push lambda amount of f_st flow along the path
            int flow_sign = +1;
            if (u < v) {
                // forward direction
                int e = graph.getEdgeId(u, v);
                table.addFlow(e, s, lambda);
                table.addFlow(e, t, -lambda);
            } else {
                // backward direction
                // push positive flow along the anti-edge
                flow_sign *= -1;

                int anti_e = graph.getEdgeId(v, u);
                table.addFlow(anti_e, s, flow_sign*lambda);
                table.addFlow(anti_e, t, -(flow_sign*lambda));
            }
        }
    }

    std::set<int> collectSubtreeVertices(std::shared_ptr<ITreeNode> node) {
        std::set<int> result(node->getMembers().begin(), node->getMembers().end());
        for (const auto& child : node->getChildren()) {
            auto sub = collectSubtreeVertices(child);
            result.insert(sub.begin(), sub.end());
        }
        return result;
    }


    void removeCycles(LinearRoutingTable& table) {
        std::set<std::pair<int, int>> all;
        for (int e = 0; e<table.src_ids.size(); e++) {
            const auto& ids  = table.src_ids[e];
            for (const auto& id : ids) {
                all.insert({id, root});
            }
        }

        for (const auto& d:all) {
            while (true) {
                auto cycle = findCycle(d, table);
                if (cycle.empty()) break;

                double minF = std::numeric_limits<double>::infinity();
                for (auto e : cycle) {
                    int e_id = graph.getEdgeId(e.first, e.second);
                    double frac = table.getFlow(e_id, d.first);
                    minF = std::min(minF,frac);
                }

                for (auto e : cycle) {
                    int e_id = graph.getEdgeId(e.first, e.second);
                    auto& dmap = table.src_ids[e_id];
                    // find index of (d.first, d.second) in table.adj_ids[e_id]
                    const auto& ids = table.src_ids[e_id];
                    size_t len = ids.size();
                    size_t lo = 0, hi = len;
                    while (lo < hi) {
                        const size_t mid = (lo + hi) >> 1;
                        const auto& mid_val = ids[mid];
                        if (mid_val < d.first)
                            lo = mid + 1;
                        else
                            hi = mid;
                    }

                    // if found, decrease the fraction
                    if (lo < len && ids[lo] == d.first) {
                        dmap[lo] -= minF;
                        if (dmap[lo] <= 0) {
                            // remove the entry
                            dmap.erase(dmap.begin() + static_cast<long>(lo));
                            table.src_flows[e_id].erase(table.src_flows[e_id].begin() + static_cast<long>(lo));
                        }
                    }

                    // also erase for the anti edge
                    int anti_e_id = graph.getEdgeId(e.second, e.first);
                    auto& rmap = table.src_ids[anti_e_id];
                    const auto& rids = table.src_ids[anti_e_id];
                    size_t rlen = rids.size();
                    size_t rlo = 0, rhi = rlen;
                    while (rlo < rhi) {
                        const size_t mid = (rlo + rhi) >> 1;
                        const auto& mid_val = rids[mid];
                        if (mid_val < d.second)
                            rlo = mid + 1;
                        else
                            rhi = mid;
                    }

                    // if found, decrease the fraction
                    if (rlo < rlen && rids[rlo] == d.second) {
                        rmap[rlo] -= minF;
                        if (rmap[rlo] <= 0) {
                            // remove the entry
                            rmap.erase(rmap.begin() + static_cast<long>(rlo));
                            table.src_flows[anti_e_id].erase(table.src_flows[anti_e_id].begin() + static_cast<long>(rlo));
                        }
                    }
                }

            }
        }
    }

    std::vector<std::pair<int, int>> findCycle(const std::pair<int, int>& d, LinearRoutingTable& table) {
        std::set<int> analyzed;
        for (int v : graph.getVertices()) {
            if (analyzed.count(v)) continue;
            std::vector<int> stack;
            if (auto maybe = findCycleRec(v, analyzed, stack, d, table))
                return *maybe;
        }
        return {};
    }

    std::optional<std::vector<std::pair<int, int>>> findCycleRec(
        int u,
        std::set<int>& analyzed,
        std::vector<int>& stack,
        const std::pair<int, int>& d,
        LinearRoutingTable& table) {
        if (std::find(stack.begin(), stack.end(), u) != stack.end())
            return std::vector<std::pair<int, int>>{};
        if (analyzed.count(u))
            return std::nullopt;

        analyzed.insert(u);
        stack.push_back(u);

        for (int w : graph.neighbors(u)) {
            int e_id = graph.getEdgeId(u, w);
            double frac = table.getFlow(e_id, d.first);
            if (frac <= 0) continue;

            if (auto child = findCycleRec(w, analyzed, stack, d, table)) {
                auto cycle = *child;
                cycle.insert(cycle.begin(), {u, w});
                return cycle;
            }
        }

        stack.pop_back();
        return std::nullopt;
    }


};



#endif //OBLIVIOUSROUTING_RAECKE_LINEAR_OR_H