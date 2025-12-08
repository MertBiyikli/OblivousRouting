//
// Created by Mert Biyikli on 18.09.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_RANDOM_MST_H
#define OBLIVIOUSROUTING_RAECKE_RANDOM_MST_H

/**
 * This class implements Räcke’s oblivious‐routing solver
 * as in his STOC 2008 paper using the uniform MST algorithm.
 */

#include <memory>
#include <vector>
#include <unordered_map>
#include <random>
#include "mst.h"
#include "../../solver/solver.h"
#include "../raecke_base.h"

// ---------------------- raecke_random_mst.h ----------------------
class RaeckeMST : public RaeckeBase<MSTTree> {
    RandomMST mst;
    std::random_device rd;
    std::mt19937_64 rng{rd()};
    std::uniform_int_distribution<uint64_t> dist;

public:
    RaeckeMST() = default;
    virtual ~RaeckeMST() = default;
    void init(GraphADJ& g);
    MSTTree getTree(GraphADJ& g);
    void computeRLoads(int idx, MSTTree& t, GraphADJ& g);

    // new version
    template <typename Transform>
    void run(Transform& transform) {
        init(m_graph);
        m_lambdaSum = 0.0;
        int id = 0;
        while (m_lambdaSum < 1.0) {
            auto start = std::chrono::high_resolution_clock::now();
            m_lambdaSum += iterate(id, transform);
            oracle_running_times.push_back(
                std::chrono::duration<double, std::milli>(
                    std::chrono::high_resolution_clock::now() - start
                ).count()
            );
            id++;
        }
    }

    template <typename Transform>
    double iterate(int idx, Transform& transform) {
        MSTTree t = getTree(m_graph);
        computeRLoads(idx, t, m_graph);
        double l = getMaxRload(idx);
        double lambda = std::min(1.0/l, 1.0 - m_lambdaSum);

        // directly add flow contribution
        transform.addTree(t, lambda, m_graph);

        // update weights
        m_lambdaSum += lambda;
        computeNewDistances(m_graph);
        return lambda;
    }

    int getIterationCount() const { return (int)oracle_running_times.size(); }
};


#endif //OBLIVIOUSROUTING_RAECKE_RANDOM_MST_H