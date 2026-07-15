#pragma once

#include <optional>
#include <unordered_map>
#include <vector>


struct HierarchyCluster {
    // Globally unique across the entire hierarchy.
    int id = -1;
    int xcut_partition_id = -1;
    // 0 = finest/leaf level, increasing toward root.
    int level = -1;

    std::vector<int> original_vertices;
    std::vector<int> induced_edges;

    std::optional<int> parent;
    std::vector<int> children;

    // Selected original graph vertex through which traffic moves upward.
    int portal = -1;
};

struct HierarchyLevel {
    // 0 = finest level.
    int level = -1;

    std::vector<HierarchyCluster> clusters;
};

struct HierarchyResult {
    std::vector<HierarchyLevel> levels;

    // Cluster ID -> level index and position inside that level.
    std::unordered_map<int, std::pair<int, int>> cluster_location;

    // Original vertex -> finest cluster containing that vertex.
    std::vector<int> vertex_to_leaf;

    std::optional<int> root;

    bool empty() const noexcept {
        return levels.empty();
    }

    const HierarchyCluster* findCluster(int id) const {
        const auto it = cluster_location.find(id);

        if (it == cluster_location.end()) {
            return nullptr;
        }

        const auto [level_index, cluster_index] = it->second;
        return &levels[level_index].clusters[cluster_index];
    }

    HierarchyCluster* findCluster(int id) {
        const auto it = cluster_location.find(id);

        if (it == cluster_location.end()) {
            return nullptr;
        }

        const auto [level_index, cluster_index] = it->second;
        return &levels[level_index].clusters[cluster_index];
    }
};