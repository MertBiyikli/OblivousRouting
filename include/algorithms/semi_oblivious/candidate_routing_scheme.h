//
// Created by Mert Biyikli on 12.06.26.
//

#ifndef OBLIVIOUSROUTING_CANDIDATE_ROUTING_SCHEME_H
#define OBLIVIOUSROUTING_CANDIDATE_ROUTING_SCHEME_H
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <functional>

struct Path {
    int source;
    int target;
    std::vector<int> edges;
    double weight = 1.0;
};

struct PairKey {
    int source;
    int target;

    bool operator==(const PairKey& other) const {
        return source == other.source && target == other.target;
    }
};

struct PairKeyHash {
    std::size_t operator()(const PairKey& key) const {
        return std::hash<int>{}(key.source) ^ (std::hash<int>{}(key.target) << 1);
    }
};

struct CandidatePathDebugStats {
    std::size_t num_pairs = 0;
    std::size_t total_paths = 0;

    double avg_paths_per_pair = 0.0;

    double avg_max_path_weight = 0.0;
    double avg_top3_weight = 0.0;
    double avg_entropy = 0.0;

    double avg_path_length = 0.0;

    double avg_pairwise_overlap = 0.0;

    double avg_unique_edges_per_pair = 0.0;
};

static double computeJaccardOverlap(
    const Path& a,
    const Path& b
) {
    std::unordered_set<int> A(
        a.edges.begin(),
        a.edges.end()
    );

    std::unordered_set<int> B(
        b.edges.begin(),
        b.edges.end()
    );

    std::size_t intersection = 0;

    for (int e : A) {
        if (B.contains(e)) {
            ++intersection;
        }
    }

    const std::size_t unionSize =
        A.size() + B.size() - intersection;

    if (unionSize == 0) {
        return 0.0;
    }

    return static_cast<double>(intersection)
         / static_cast<double>(unionSize);
}

class CandidateRoutingScheme {
public:
    void addPath(int s, int t, Path path) {
        paths_[{s, t}].push_back(std::move(path));
    }

    [[nodiscard]] bool hasPaths(int s, int t) const {
        return paths_.contains({s, t});
    }

    [[nodiscard]] const std::vector<Path>& paths(int s, int t) const {
        auto it = paths_.find({s, t});
        if (it == paths_.end()) {
            throw std::runtime_error("No candidate paths for demand pair");
        }
        return it->second;
    }

    [[nodiscard]] const auto& allPaths() const {
        return paths_;
    }

    [[nodiscard]] std::size_t numPaths() const {
        std::size_t total = 0;

        for (const auto& [pair, paths] : paths_) {
            total += paths.size();
        }

        return total;
    }

    [[nodiscard]] double averagePathsPerPair(int n) const {
        if (n <= 1) {
            return 0.0;
        }

        return static_cast<double>(numPaths()) /
               static_cast<double>(n * (n - 1));
    }

    [[nodiscard]] std::size_t numPairsWithPaths() const {
        return paths_.size();
    }

    CandidatePathDebugStats computeDebugStats() const {
    CandidatePathDebugStats stats;
    stats.num_pairs = paths_.size();

    if (stats.num_pairs == 0) {
        return stats;
    }

    double sumMaxWeight = 0.0;
    double sumTop3Weight = 0.0;
    double sumEntropy = 0.0;

    double sumPathLength = 0.0;
    std::size_t totalPathsSeen = 0;

    double sumPairwiseOverlap = 0.0;
    std::size_t overlapComparisons = 0;

    double sumUniqueEdgesPerPair = 0.0;

    for (const auto& [pair, paths] : paths_) {
        stats.total_paths += paths.size();

        std::unordered_set<int> pairEdges;

        for (const auto& path : paths) {
            sumPathLength += static_cast<double>(path.edges.size());
            ++totalPathsSeen;

            for (int e : path.edges) {
                pairEdges.insert(e);
            }
        }

        sumUniqueEdgesPerPair += static_cast<double>(pairEdges.size());

        for (std::size_t i = 0; i < paths.size(); ++i) {
            for (std::size_t j = i + 1; j < paths.size(); ++j) {
                sumPairwiseOverlap += computeJaccardOverlap(paths[i], paths[j]);
                ++overlapComparisons;
            }
        }

        double totalWeight = 0.0;
        std::vector<double> weights;

        for (const auto& path : paths) {
            const double w = std::max(0.0, path.weight);
            totalWeight += w;
            weights.push_back(w);
        }

        if (totalWeight <= 1e-12) {
            continue;
        }

        std::sort(weights.begin(), weights.end(), std::greater<>());

        sumMaxWeight += weights.front() / totalWeight;

        double top3 = 0.0;
        for (std::size_t i = 0; i < std::min<std::size_t>(3, weights.size()); ++i) {
            top3 += weights[i];
        }

        sumTop3Weight += top3 / totalWeight;

        double entropy = 0.0;
        for (double w : weights) {
            const double p = w / totalWeight;
            if (p > 1e-12) {
                entropy -= p * std::log(p);
            }
        }

        sumEntropy += entropy;
    }

    stats.avg_paths_per_pair =
        static_cast<double>(stats.total_paths) /
        static_cast<double>(stats.num_pairs);

    stats.avg_max_path_weight =
        sumMaxWeight / static_cast<double>(stats.num_pairs);

    stats.avg_top3_weight =
        sumTop3Weight / static_cast<double>(stats.num_pairs);

    stats.avg_entropy =
        sumEntropy / static_cast<double>(stats.num_pairs);

    if (totalPathsSeen > 0) {
        stats.avg_path_length =
            sumPathLength / static_cast<double>(totalPathsSeen);
    }

    if (overlapComparisons > 0) {
        stats.avg_pairwise_overlap =
            sumPairwiseOverlap / static_cast<double>(overlapComparisons);
    }

    stats.avg_unique_edges_per_pair =
        sumUniqueEdgesPerPair / static_cast<double>(stats.num_pairs);

    return stats;
}

private:
    std::unordered_map<PairKey, std::vector<Path>, PairKeyHash> paths_;
};

inline void printCandidatePathDebugStats(
    const std::string& routingBase,
    const CandidateRoutingScheme& scheme
) {
    const auto stats = scheme.computeDebugStats();

    std::cout << "Candidate path diagnostics [" << routingBase << "]\n";
    std::cout << "  Pairs with paths: " << stats.num_pairs << '\n';
    std::cout << "  Total paths: " << stats.total_paths << '\n';
    std::cout << "  Avg paths/pair: " << stats.avg_paths_per_pair << '\n';
    std::cout << "  Avg max path weight: " << stats.avg_max_path_weight << '\n';
    std::cout << "  Avg top-3 path weight: " << stats.avg_top3_weight << '\n';
    std::cout << "  Avg entropy: " << stats.avg_entropy << '\n';
    std::cout
    << "  Avg path length: "
    << stats.avg_path_length
    << '\n';

    std::cout
        << "  Avg pairwise overlap: "
        << stats.avg_pairwise_overlap
        << '\n';

    std::cout
        << "  Avg unique edges/pair: "
        << stats.avg_unique_edges_per_pair
        << '\n';
}

#endif //OBLIVIOUSROUTING_CANDIDATE_ROUTING_SCHEME_H