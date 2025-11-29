//
// Created by Mert Biyikli on 25.10.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_CKR_TRANSFORM_H
#define OBLIVIOUSROUTING_RAECKE_CKR_TRANSFORM_H


#include <queue>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include "../../utils/hash.h"
#include "../raecke_transform_base.h"
#include "ckr_tree_decomposer.h"
#include "../../datastructures/graph_csr.h"
#include "../raecke_oracle_iteration.h"

class RaeckeCKRTransform : public RaeckeTransformBase<std::shared_ptr<TreeNode>> {
public:
    virtual EdgeDemandMap& addTree(std::shared_ptr<TreeNode>& root, double lambda, Graph& graph) override{
        distributeDemands(root, graph, lambda);
        normalizeOldSolutionBasedOnNewLambda(lambda);
        removeCycles(graph);
        return arc2demand2cumulativeFraction;
    }

    virtual EdgeDemandMap& addTree(std::shared_ptr<TreeNode>& root, double lambda, Graph_csr& graph) {
        distributeDemands(root, graph, lambda);
        normalizeOldSolutionBasedOnNewLambda(lambda);
        removeCycles(graph);
        return arc2demand2cumulativeFraction;
    }



private:
    void distributeDemands(const std::shared_ptr<TreeNode>& node, Graph& graph, double lambda) {
        if (!node) return;

        for (const auto& child : node->children) {
            std::set<int> A = collectSubtreeVertices(child);
            std::set<int> B;
            for (int v : graph.getVertices()) {
                if (!A.count(v)) B.insert(v);
            }

            for (int src : A) {
                for (int dst : B) {
                    TerminalPair pair{src, dst};
                    TerminalPair revPair{dst, src};
                    auto path = graph.getShortestPath(src, dst);

                    for (size_t i = 0; i + 1 < path.size(); ++i) {
                        std::pair<int, int> arc{path[i], path[i+1]};

                        auto& fwdMap = arc2demand2cumulativeFraction[arc];
                        auto& revMap = arc2demand2cumulativeFraction[{arc.second, arc.first}];

                        fwdMap[pair] += lambda;
                        revMap[revPair] += lambda;
                    }
                }
            }

            distributeDemands(child, graph, lambda);
        }
    }

    void distributeDemands(const std::shared_ptr<TreeNode>& node, Graph_csr& graph, double lambda) {
        if (!node) return;

        for (const auto& child : node->children) {
            std::set<int> A = collectSubtreeVertices(child);
            std::set<int> B;
            for (int v : graph.vertices) {
                if (!A.count(v)) B.insert(v);
            }

            for (int src : A) {
                for (int dst : B) {
                    TerminalPair pair{src, dst};
                    TerminalPair revPair{dst, src};
                    auto path = graph.getShortestPath(src, dst);

                    for (size_t i = 0; i + 1 < path.size(); ++i) {
                        std::pair<int, int> arc{path[i], path[i+1]};

                        auto& fwdMap = arc2demand2cumulativeFraction[arc];
                        auto& revMap = arc2demand2cumulativeFraction[{arc.second, arc.first}];

                        fwdMap[pair] += lambda;
                        revMap[revPair] += lambda;
                    }
                }
            }

            distributeDemands(child, graph, lambda);
        }
    }

    std::set<int> collectSubtreeVertices(const std::shared_ptr<TreeNode>& node) {
        std::set<int> result(node->members.begin(), node->members.end());
        for (const auto& child : node->children) {
            auto sub = collectSubtreeVertices(child);
            result.insert(sub.begin(), sub.end());
        }
        return result;
    }

    void normalizeOldSolutionBasedOnNewLambda(double lambda) override{
        for (auto& [arc, dmap] : arc2demand2cumulativeFraction) {
            for (auto& [d, frac] : dmap) {
                frac *= (1.0 - lambda);
            }
        }
    }

    void removeCycles(Graph& graph) {
        std::set<std::pair<int, int>> allD;
        for (auto& [e, dmap] : arc2demand2cumulativeFraction)
            for (auto& [d, _] : dmap)
                allD.insert({d.first, d.second});

        for (auto d : allD) {
            while (true) {
                auto cycle = findCycle(d, graph);
                if (cycle.empty()) break;

                double minF = std::numeric_limits<double>::infinity();
                for (auto e : cycle)
                    minF = std::min(minF, arc2demand2cumulativeFraction[e][d]);

                for (auto e : cycle) {
                    auto& dmap = arc2demand2cumulativeFraction[e];
                    dmap[d] -= minF;
                    if (dmap[d] <= 0) dmap.erase(d);
                    if (dmap.empty()) arc2demand2cumulativeFraction.erase(e);

                    auto& rmap = arc2demand2cumulativeFraction[{e.second, e.first}];
                    rmap[{d.second, d.first}] -= minF;
                    if (rmap[{d.second, d.first}] <= 0) rmap.erase({d.second, d.first});
                    if (rmap.empty()) arc2demand2cumulativeFraction.erase({e.second, e.first});
                }
            }
        }
    }

    std::vector<std::pair<int, int>> findCycle(const std::pair<int, int>& d, Graph& graph) {
        std::set<int> analyzed;
        for (int v : graph.getVertices()) {
            if (analyzed.count(v)) continue;
            std::vector<int> stack;
            if (auto maybe = findCycleRec(v, analyzed, stack, d, graph))
                return *maybe;
        }
        return {};
    }

    std::optional<std::vector<std::pair<int, int>>> findCycleRec(
        int u,
        std::set<int>& analyzed,
        std::vector<int>& stack,
        const std::pair<int, int>& d,
        Graph& graph) {

        if (std::find(stack.begin(), stack.end(), u) != stack.end())
            return std::vector<std::pair<int, int>>{};
        if (analyzed.count(u))
            return std::nullopt;

        analyzed.insert(u);
        stack.push_back(u);

        for (int w : graph.neighbors(u)) {
            std::pair<int, int> e{u, w};
            auto it = arc2demand2cumulativeFraction.find(e);
            if (it == arc2demand2cumulativeFraction.end()) continue;
            if (!it->second.count(d)) continue;

            if (auto child = findCycleRec(w, analyzed, stack, d, graph)) {
                auto cycle = *child;
                cycle.insert(cycle.begin(), e);
                return cycle;
            }
        }

        stack.pop_back();
        return std::nullopt;
    }

     void removeCycles(Graph_csr& graph) {
        std::set<std::pair<int, int>> allD;
        for (auto& [e, dmap] : arc2demand2cumulativeFraction)
            for (auto& [d, _] : dmap)
                allD.insert({d.first, d.second});

        for (auto d : allD) {
            while (true) {
                auto cycle = findCycle(d, graph);
                if (cycle.empty()) break;

                double minF = std::numeric_limits<double>::infinity();
                for (auto e : cycle)
                    minF = std::min(minF, arc2demand2cumulativeFraction[e][d]);

                for (auto e : cycle) {
                    auto& dmap = arc2demand2cumulativeFraction[e];
                    dmap[d] -= minF;
                    if (dmap[d] <= 0) dmap.erase(d);
                    if (dmap.empty()) arc2demand2cumulativeFraction.erase(e);

                    auto& rmap = arc2demand2cumulativeFraction[{e.second, e.first}];
                    rmap[{d.second, d.first}] -= minF;
                    if (rmap[{d.second, d.first}] <= 0) rmap.erase({d.second, d.first});
                    if (rmap.empty()) arc2demand2cumulativeFraction.erase({e.second, e.first});
                }
            }
        }
    }

    std::vector<std::pair<int, int>> findCycle(const std::pair<int, int>& d, Graph_csr& graph) {
        std::set<int> analyzed;
        for (int v : graph.vertices) {
            if (analyzed.count(v)) continue;
            std::vector<int> stack;
            if (auto maybe = findCycleRec(v, analyzed, stack, d, graph))
                return *maybe;
        }
        return {};
    }

    std::optional<std::vector<std::pair<int, int>>> findCycleRec(
        int u,
        std::set<int>& analyzed,
        std::vector<int>& stack,
        const std::pair<int, int>& d,
        Graph_csr& graph) {

        if (std::find(stack.begin(), stack.end(), u) != stack.end())
            return std::vector<std::pair<int, int>>{};
        if (analyzed.count(u))
            return std::nullopt;

        analyzed.insert(u);
        stack.push_back(u);

        for (auto w : graph.neighbors(u)) {
            std::pair<int, int> e{u, w};
            auto it = arc2demand2cumulativeFraction.find(e);
            if (it == arc2demand2cumulativeFraction.end()) continue;
            if (!it->second.count(d)) continue;

            if (auto child = findCycleRec(w, analyzed, stack, d, graph)) {
                auto cycle = *child;
                cycle.insert(cycle.begin(), e);
                return cycle;
            }
        }

        stack.pop_back();
        return std::nullopt;
    }
};




#endif //OBLIVIOUSROUTING_RAECKE_CKR_TRANSFORM_H