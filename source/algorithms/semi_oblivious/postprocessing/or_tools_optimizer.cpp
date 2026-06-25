//
// Created by Mert Biyikli on 12.06.26.
//

#include "../../../../include/algorithms/semi_oblivious/postprocessing/or_tools_optimizer.h"
#include "routing/storage/allpair_routing_table.h"
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "ortools/linear_solver/linear_solver.h"

using operations_research::MPSolver;
using operations_research::MPVariable;
using operations_research::MPConstraint;
using operations_research::MPObjective;

struct PairPathVarKey {
    int source;
    int target;
    std::size_t pathIndex;

    bool operator==(const PairPathVarKey& other) const {
        return source == other.source &&
               target == other.target &&
               pathIndex == other.pathIndex;
    }
};

struct PairPathVarKeyHash {
    std::size_t operator()(const PairPathVarKey& key) const {
        std::size_t h1 = std::hash<int>{}(key.source);
        std::size_t h2 = std::hash<int>{}(key.target);
        std::size_t h3 = std::hash<std::size_t>{}(key.pathIndex);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};
struct SemiObliviousLPDiagnostics {
    std::size_t activeDemandPairs = 0;
    std::size_t totalVariables = 0;
    std::size_t activeVariables = 0;
    std::size_t splitPairs = 0;
    std::size_t totalActivePaths = 0;
    std::size_t maxActivePathsForPair = 0;

    double avgLargestAlpha = 0.0;
    double avgTop3Alpha = 0.0;
    double avgAlphaEntropy = 0.0;

    [[nodiscard]] double activePathRatio() const {
        if (totalVariables == 0) {
            return 0.0;
        }

        return static_cast<double>(activeVariables) /
               static_cast<double>(totalVariables);
    }

    [[nodiscard]] double avgActivePathsPerPair() const {
        if (activeDemandPairs == 0) {
            return 0.0;
        }

        return static_cast<double>(totalActivePaths) /
               static_cast<double>(activeDemandPairs);
    }
};
SemiObliviousLPDiagnostics computeLPDiagnostics(
    const CandidateRoutingScheme& candidateScheme,
    const demands& demand,
    const std::unordered_map<PairPathVarKey, MPVariable*, PairPathVarKeyHash>& xVars
) {
    constexpr double activeTol = 1e-6;

    SemiObliviousLPDiagnostics diagnostics;

    double sumLargestAlpha = 0.0;
    double sumTop3Alpha = 0.0;
    double sumAlphaEntropy = 0.0;

    for (const auto& [pair, paths] : candidateScheme.allPaths()) {
        const int s = pair.source;
        const int t = pair.target;

        const auto dOpt = demand.getDemandValue(s, t);

        if (!dOpt || *dOpt <= 0.0) {
            continue;
        }

        std::vector<double> alphas;
        std::size_t activeForPair = 0;

        for (std::size_t p = 0; p < paths.size(); ++p) {
            auto it = xVars.find(PairPathVarKey{s, t, p});

            if (it == xVars.end()) {
                continue;
            }

            const double alpha = it->second->solution_value();

            alphas.push_back(alpha);
            ++diagnostics.totalVariables;

            if (alpha > activeTol) {
                ++diagnostics.activeVariables;
                ++activeForPair;
            }
        }

        if (alphas.empty()) {
            continue;
        }

        ++diagnostics.activeDemandPairs;

        std::sort(alphas.begin(), alphas.end(), std::greater<>());

        sumLargestAlpha += alphas.front();

        double top3 = 0.0;
        for (std::size_t i = 0; i < std::min<std::size_t>(3, alphas.size()); ++i) {
            top3 += alphas[i];
        }

        sumTop3Alpha += top3;

        double entropy = 0.0;
        for (double alpha : alphas) {
            if (alpha > 1e-12) {
                entropy -= alpha * std::log(alpha);
            }
        }

        sumAlphaEntropy += entropy;

        if (activeForPair > 1) {
            ++diagnostics.splitPairs;
        }

        diagnostics.totalActivePaths += activeForPair;
        diagnostics.maxActivePathsForPair =
            std::max(diagnostics.maxActivePathsForPair, activeForPair);
    }

    if (diagnostics.activeDemandPairs > 0) {
        const double denom =
            static_cast<double>(diagnostics.activeDemandPairs);

        diagnostics.avgLargestAlpha = sumLargestAlpha / denom;
        diagnostics.avgTop3Alpha = sumTop3Alpha / denom;
        diagnostics.avgAlphaEntropy = sumAlphaEntropy / denom;
    }

    return diagnostics;
}
void printLPDiagnostics(
    const SemiObliviousLPDiagnostics& d,
    double lambda,
    int numVariables,
    int numConstraints
) {
    std::cout << "\n";
    std::cout << "Semi-oblivious LP diagnostics\n";
    std::cout << "----------------------------------------\n";
    std::cout << "  Objective lambda:              " << lambda << '\n';
    std::cout << "  LP variables:                  " << numVariables << '\n';
    std::cout << "  LP constraints:                " << numConstraints << '\n';
    std::cout << "  Active demand pairs:           " << d.activeDemandPairs << '\n';
    std::cout << "  Active path variables:         "
              << d.activeVariables << " / " << d.totalVariables << '\n';
    std::cout << "  Active path ratio:             "
              << d.activePathRatio() << '\n';
    std::cout << "  Split pairs:                   "
              << d.splitPairs << " / " << d.activeDemandPairs << '\n';
    std::cout << "  Avg largest alpha:             "
              << d.avgLargestAlpha << '\n';
    std::cout << "  Avg top-3 alpha mass:          "
              << d.avgTop3Alpha << '\n';
    std::cout << "  Avg alpha entropy:             "
              << d.avgAlphaEntropy << '\n';
    std::cout << "  Avg active paths / pair:       "
              << d.avgActivePathsPerPair() << '\n';
    std::cout << "  Max active paths for one pair: "
              << d.maxActivePathsForPair << '\n';
    std::cout << "----------------------------------------\n";
}
void validateSemiObliviousExtraction(
    const IGraph& graph,
    const CandidateRoutingScheme& candidateScheme,
    const demands& demand,
    const std::unordered_map<PairPathVarKey, MPVariable*, PairPathVarKeyHash>& xVars,
    const RoutingScheme& extractedScheme,
    double lambda
) {
    const int m = graph.getNumDirectedEdges();
    const double tol = 1e-6;

    std::vector<double> lpLoad(m, 0.0);
    std::vector<double> schemeLoad(m, 0.0);

    // 1. Compute LP-implied load.
    for (const auto& [pair, paths] : candidateScheme.allPaths()) {
        const int s = pair.source;
        const int t = pair.target;

        const auto dOpt = demand.getDemandValue(s, t);
        if (!dOpt || *dOpt <= 0.0) {
            continue;
        }

        const double d = *dOpt;

        for (std::size_t p = 0; p < paths.size(); ++p) {
            auto it = xVars.find(PairPathVarKey{s, t, p});
            if (it == xVars.end()) {
                continue;
            }

            const double alpha = it->second->solution_value();
            if (std::abs(alpha) <= tol) {
                continue;
            }

            for (int e : paths[p].edges) {
                lpLoad[e] += d * alpha;
            }
        }
    }

    // 2. Compute extracted-scheme load.
    for (int s = 0; s < graph.getNumNodes(); ++s) {
        for (int t = 0; t < graph.getNumNodes(); ++t) {
            if (s == t) {
                continue;
            }

            const auto dOpt = demand.getDemandValue(s, t);
            if (!dOpt || *dOpt <= 0.0) {
                continue;
            }

            const double d = *dOpt;

            for (int e = 0; e < m; ++e) {
                const double f = extractedScheme.getFlow(e, s, t);
                schemeLoad[e] += d * std::abs(f);
            }
        }
    }

    // 3. Compare canonical undirected loads.
    double maxDiff = 0.0;
    int worstEdge = -1;

    double maxLpCongestion = 0.0;
    double maxSchemeCongestion = 0.0;

    std::vector<bool> processed(m, false);

    for (int e = 0; e < m; ++e) {
        if (processed[e]) {
            continue;
        }

        const int anti = graph.getAntiEdge(e);

        processed[e] = true;
        if (anti >= 0) {
            processed[anti] = true;
        }

        const double lpUndirectedLoad =
            lpLoad[e] + (anti >= 0 ? lpLoad[anti] : 0.0);

        const double schemeUndirectedLoad =
            schemeLoad[e] + (anti >= 0 ? schemeLoad[anti] : 0.0);

        const double diff = std::abs(lpUndirectedLoad - schemeUndirectedLoad);

        if (diff > maxDiff) {
            maxDiff = diff;
            worstEdge = e;
        }

        const double capacity = graph.getEdgeCapacity(e);

        if (capacity > 0.0) {
            maxLpCongestion =
                std::max(maxLpCongestion, lpUndirectedLoad / capacity);

            maxSchemeCongestion =
                std::max(maxSchemeCongestion, schemeUndirectedLoad / capacity);
        }
    }

    std::cout << "[SemiObliviousValidation] lambda = "
              << lambda << '\n';

    std::cout << "[SemiObliviousValidation] LP congestion = "
              << maxLpCongestion << '\n';

    std::cout << "[SemiObliviousValidation] extracted scheme congestion = "
              << maxSchemeCongestion << '\n';

    std::cout << "[SemiObliviousValidation] max load diff = "
              << maxDiff
              << " on edge "
              << worstEdge
              << '\n';
}

SemiObliviousOptimizationResult OrToolsSemiObliviousLoadOptimizer::optimize(
    const IGraph& graph,
    const CandidateRoutingScheme& candidateScheme,
    const demands& _demand
) {
    MPSolver solver(
        "semi_oblivious_load_optimizer",
        MPSolver::GLOP_LINEAR_PROGRAMMING
    );

    const double inf = solver.infinity();

    auto* lambda = solver.MakeNumVar(0.0, inf, "lambda");

    std::unordered_map<PairPathVarKey, MPVariable*, PairPathVarKeyHash> xVars;

    // 1. Create path split variables and demand split constraints.
    for (const auto& [pair, paths] : candidateScheme.allPaths()) {
        const int s = pair.source;
        const int t = pair.target;

        auto d = _demand.getDemandValue(s, t); // adapt if your API differs

        if (!d || *d <= 0.0) {
            continue;
        }

        if (paths.empty()) {
            throw std::runtime_error("Demand pair has no candidate paths");
        }

        auto* splitConstraint = solver.MakeRowConstraint(1.0, 1.0);

        for (std::size_t p = 0; p < paths.size(); ++p) {
            const std::string name =
                "x_" + std::to_string(s) + "_" +
                std::to_string(t) + "_" +
                std::to_string(p);

            auto* x = solver.MakeNumVar(0.0, inf, name);

            xVars.emplace(PairPathVarKey{s, t, p}, x);
            splitConstraint->SetCoefficient(x, 1.0);
        }
    }

    // 2. Edge congestion constraints:
    //
    // sum_{s,t,p: e in p} d_st * x_stp <= lambda * capacity(e)
    // 2. Undirected edge congestion constraints:
    //
    // load(e) + load(reverse(e)) <= lambda * capacity(e)
    //
    // This must match computeRoutingSchemeCongestion if the evaluator aggregates
    // both directions of an undirected edge.
    std::vector<bool> processed(graph.getNumDirectedEdges(), false);

    for (int e = 0; e < graph.getNumDirectedEdges(); ++e) {
        if (processed[e]) {
            continue;
        }

        const int rev = graph.getAntiEdge(e); // adapt name if needed

        processed[e] = true;
        if (rev >= 0) {
            processed[rev] = true;
        }

        const double capacity = graph.getEdgeCapacity(e);

        if (capacity <= 0.0) {
            throw std::runtime_error("Edge has non-positive capacity");
        }

        auto* edgeConstraint = solver.MakeRowConstraint(-inf, 0.0);
        edgeConstraint->SetCoefficient(lambda, -capacity);

        for (const auto& [pair, paths] : candidateScheme.allPaths()) {
            const int s = pair.source;
            const int t = pair.target;

            const auto d = _demand.getDemandValue(s, t);

            if (!d || *d <= 0.0) {
                continue;
            }

            for (std::size_t p = 0; p < paths.size(); ++p) {
                const auto& path = paths[p];

                bool usesUndirectedEdge = false;

                for (int pathEdge : path.edges) {
                    if (pathEdge == e || pathEdge == rev) {
                        usesUndirectedEdge = true;
                        break;
                    }
                }

                if (!usesUndirectedEdge) {
                    continue;
                }

                auto it = xVars.find(PairPathVarKey{s, t, p});
                if (it == xVars.end()) {
                    throw std::runtime_error("Missing LP variable for path");
                }

                edgeConstraint->SetCoefficient(it->second, *d);
            }
        }
    }

    // 3. Objective: minimize lambda.
    MPObjective* objective = solver.MutableObjective();
    objective->SetCoefficient(lambda, 1.0);
    objective->SetMinimization();

    const auto status = solver.Solve();

    if (status != MPSolver::OPTIMAL
        && status != MPSolver::FEASIBLE) {
        throw std::runtime_error("Semi-oblivious LP infeasible or not solved");
    }

    const auto diagnostics = computeLPDiagnostics(candidateScheme,_demand,xVars);


    AllPairRoutingTable table;
    table.init(graph);

    for (const auto& [pair, paths] : candidateScheme.allPaths()) {
        const int s = pair.source;
        const int t = pair.target;

        const auto d = _demand.getDemandValue(s, t);

        if (!d || *d <= 0.0) {
            continue;
        }

        for (std::size_t p = 0; p < paths.size(); ++p) {
            auto it = xVars.find(PairPathVarKey{s, t, p});

            if (it == xVars.end()) {
                throw std::runtime_error("Missing LP variable during solution extraction");
            }

            const double alpha = it->second->solution_value();

            if (std::abs(alpha) <= 1e-9) {
                continue;
            }

            for (int e : paths[p].edges) {
                table.addFlow(e, s, t, alpha);
            }
        }
    }

    SemiObliviousOptimizationResult result;
     result.scheme = std::make_unique<AllPairRoutingScheme>(graph, std::move(table));


    return result;
}