//
// Created by Mert Biyikli on 12.06.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_TRANSFORM_H
#define OBLIVIOUSROUTING_RAECKE_TRANSFORM_H

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include "../utils/hash.h"
#include "raecke_tree_decomp.h"
// RaeckeSolutionTransform.hpp
#pragma once

#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <utility>
#include "../utils/hash.h"

class RaeckeSolutionTransform {
public:
    // Recorded decomposition trees (if needed)
    std::vector<std::set<std::pair<int, int> >> decTree2notSpanningTrees;
    // Routing: arcID -> ( (src,dst) -> cumulative fraction )
    EdgeDemandMap arc2demand2cumulativeFraction;

    RaeckeSolutionTransform() = default;

    const EdgeDemandMap& getRoutingTable() const {
        return arc2demand2cumulativeFraction;
    }

    /**
     * Add a per-destination set of trees (vertex-disjoint) with weight lambda.
     * destination2source2path: destID -> ( srcID -> set<ArcID> )
     * Returns the updated arc2demand2cumulativeFraction.
     */
    EdgeDemandMap& addASetOfTreesPerDestination(
            const std::map<int, std::map<int, std::set<std::pair<int, int>>>>& destination2source2path,
            double lambda)
    {
        normalizeOldSolutionBasedOnNewLambda(lambda);

        std::cout << "Adding a set of trees with lambda = " << lambda << "\n";
        for (auto const& [dest, srcMap] : destination2source2path) {
            for (auto const& [src, arcSet] : srcMap) {
                if(dest > src) continue; // skip trivial cases
                std::cout << "Adding trees for destination " << dest << " from source " << src << "\n";

                TerminalPair demand{src, dest};
                for (std::pair<int, int> arcId : arcSet) {
                    auto& demand2frac = arc2demand2cumulativeFraction[arcId];
                    demand2frac[demand] += lambda;
                }
            }
        }
        return arc2demand2cumulativeFraction;
    }

    /**
     * Add one FRT decomposition tree with weight lambda, updating demands.
     * Mirrors Java addTree(Tree t, Double lambda, Graph graph).
     */
    EdgeDemandMap& addTree(
            FRT_Tree& t,
            double           lambda,
            Graph&    graph)
    {
        normalizeOldSolutionBasedOnNewLambda(lambda);

        // Record this tree's arcs (optional)
        std::set<std::pair<int, int> > tree;
        decTree2notSpanningTrees.push_back(tree);
        auto& treeArcs = decTree2notSpanningTrees.back();

        // BFS on decomposition tree
        std::queue<std::shared_ptr<FRT_Node>> q;
        for (auto& c : t.GetRoot()->getChildren())
            q.push(c);

        while (!q.empty()) {
            auto node   = q.front(); q.pop();
            auto parent = node->getParent();

            // skip trivial splits
            if (parent->getCenter() == node->getCenter()) {
                for (auto& c : node->getChildren())
                    q.push(c);
                continue;
            }

            // Partition vertices: A = node's vertices, B = complement
            auto A = node->GetVertices();
            auto B = graph.getVertices();
            for (int v : A) {
                auto it = std::find(B.begin(), B.end(), v);
                if (it != B.end()) B.erase(it);
            }

            // For each demand (src in A, dst in B)
            for (int src : A) {
                for (int dst : B) {
                    TerminalPair pair{src, dst};
                    TerminalPair revPair{dst, src};

                    // Walk the SP from src‐center to dst‐center
                    auto path = graph.getShortestPath(
                            node->getCenter(),
                            parent->getCenter());

                    for (size_t i = 0; i + 1 < path.size(); ++i) {


                        // record arc in this tree
                        treeArcs.insert({path[i], path[i+1]});

                        // update directed demands
                        auto& fwdMap = arc2demand2cumulativeFraction[{path[i], path[i+1]}];
                        auto& revMap = arc2demand2cumulativeFraction[{path[i+1], path[i]}];

                        double prev = 0.0;
                        if (auto it = fwdMap.find(pair); it != fwdMap.end()) prev = it->second;

                        fwdMap[pair] = prev + lambda;
                        revMap[revPair] = prev + lambda;
                    }
                }
            }

            // enqueue children
            for (auto& c : node->getChildren())
                q.push(c);
        }

        removeCycles(graph);
        return arc2demand2cumulativeFraction;
    }

private:
    void normalizeOldSolutionBasedOnNewLambda(double lambda) {
        for (auto& [arc, dmap] : arc2demand2cumulativeFraction) {
            for (auto& [d, frac] : dmap) {
                frac *= (1.0 - lambda);
            }
        }
    }

    // Remove cycles for each directed demand
    void removeCycles(Graph& graph) {
        // gather demands
        std::set<std::pair<int, int>> allD;
        for (auto& [e, dmap] : arc2demand2cumulativeFraction)
            for (auto& [d, _] : dmap)
                allD.insert({d.first, d.second});

        for (auto d : allD) {
            while (true) {
                auto cycle = findCycle(d, graph);
                if (cycle.empty()) break;

                // find bottleneck
                double minF = std::numeric_limits<double>::infinity();
                for (auto e : cycle) {
                    minF = std::min(minF, arc2demand2cumulativeFraction[e][d]);
                }
                std::pair<int, int> rd{d.first, d.second};

                // subtract
                for (auto e : cycle) {
                    auto& dmap = arc2demand2cumulativeFraction[e];
                    dmap[{d.first, d.second}] -= minF;
                    if (dmap[{d.first, d.second}] <= 0) dmap.erase({d.first, d.second});
                    if (dmap.empty()) arc2demand2cumulativeFraction.erase(e);

                    auto& rmap = arc2demand2cumulativeFraction[{e.second, e.first}];
                    rmap[rd] -= minF;
                    if (rmap[rd] <= 0) rmap.erase(rd);
                    if (rmap.empty()) arc2demand2cumulativeFraction.erase({e.second, e.first});
                }
            }
        }
    }

    // Find one cycle carrying demand d (directed)
    std::vector<std::pair<int, int> > findCycle(
            const std::pair<int, int>&      d,
            Graph&       graph)
    {
        std::set<int> analyzed;
        for (int v : graph.getVertices()) {
            if (analyzed.count(v)) continue;
            std::vector<int> stack;
            if (auto maybe = findCycleRec(v, analyzed, stack, d, graph))
                return *maybe;
        }
        return {};
    }

    // DFS helper
    std::optional<std::vector<std::pair<int, int>>> findCycleRec(
            int                  u,
            std::set<int>&       analyzed,
            std::vector<int>&    stack,
            const std::pair<int, int>&        d,
            Graph&         graph)
    {
        if (std::find(stack.begin(), stack.end(), u) != stack.end())
            return std::vector<std::pair<int, int>>{};
        if (analyzed.count(u))
            return std::nullopt;

        analyzed.insert(u);
        stack.push_back(u);

        for (int w : graph.neighbors(u)) {
            std::pair<int, int> e={u, w};
            auto it = arc2demand2cumulativeFraction.find(e);
            if (it == arc2demand2cumulativeFraction.end()) continue;
            auto& dmap = it->second;
            if (!dmap.count(d)) continue;

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

#endif //OBLIVIOUSROUTING_RAECKE_TRANSFORM_H
