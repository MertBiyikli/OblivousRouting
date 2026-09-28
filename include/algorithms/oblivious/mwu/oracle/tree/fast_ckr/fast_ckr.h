#ifndef OBLIVIOUSROUTING_FAST_CKR_H
#define OBLIVIOUSROUTING_FAST_CKR_H

#include "../tree_oracle.h"

#include <cassert>
#include <functional>
#include <limits>
#include <queue>
#include <random>
#include <utility>
#include <vector>

template <typename T>
class FastCKR : public TreeOracle<T>
{
public:
    explicit FastCKR(optimized::Graph<EdgeData>& g)
        : TreeOracle<T>(g) {
    }

    explicit FastCKR(optimized::Graph<EdgeData>& g, bool mendelscaling)
        : TreeOracle<T>(g, mendelscaling) {
    }

    void computeLevelPartition(optimized::Graph<EdgeData>& _g,HSTLevel& level,const std::vector<int>& x_perm,double delta) override {
        const int n = _g.getNumNodes();

        // Sample R in [Δ/4, Δ/2]
        std::mt19937_64 gen(std::random_device{}());
        std::uniform_real_distribution<double> dist(delta / 4.0, delta / 2.0);
        const double R = dist(gen);

        level.R = R;
        level.owner.assign(n, -1);
        level.centers.clear();

        std::vector<double> estimated_distances(
            n,
            std::numeric_limits<double>::infinity()
        );

        std::vector<int> P(n, 0);

        // Min-priority queue: (distance, node)
        using QueueEntry = std::pair<double, int>;

        std::priority_queue<
            QueueEntry,
            std::vector<QueueEntry>,
            std::greater<QueueEntry>
        > Q;

        for (std::size_t i = 0; i < x_perm.size(); ++i)
        {
            const int source = x_perm[i];

            if (level.owner[source] != -1)
            {
                continue;
            }

            // Unassigned vertex becomes a new center.
            level.centers.push_back(source);
            P[source] = -1;

            double& dist_source = estimated_distances[source];

            if (dist_source > 0.0)
            {
                dist_source = 0.0;
                Q.emplace(0.0, source);
            }

            while (!Q.empty() && Q.top().first <= R)
            {
                const auto [dist_w, w] = Q.top();
                Q.pop();

                // std::priority_queue has no decrease-key.
                // Ignore entries superseded by a shorter distance.
                if (dist_w > estimated_distances[w])
                {
                    continue;
                }

                if (level.owner[w] != -1)
                {
                    continue;
                }

                level.owner[w] = source;

                // Label assignment
                if (P[w] == 0)
                {
                    P[w] = static_cast<int>(i) + 1;
                }

                if (_g.edgesOf(w).empty())
                {
                    continue;
                }

                for (const auto& e : _g.edgesOf(w))
                {
                    const int u = e.tail;

                    if (u < 0 || u >= n)
                    {
                        continue;
                    }

                    if (level.owner[u] != -1)
                    {
                        continue;
                    }

                    const double new_dist =
                        dist_w + _g.edgeData(e.id).weight;

                    if (new_dist < estimated_distances[u])
                    {
                        estimated_distances[u] = new_dist;

                        // No decrease-key: insert another entry.
                        Q.emplace(new_dist, u);
                    }
                }
            }
        }

        assert(level.centers.size() <= level.owner.size());
    }
};

#endif // OBLIVIOUSROUTING_FAST_CKR_H
