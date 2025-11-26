//
// Created by Mert Biyikli on 18.09.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_MST_SOLVER_H
#define OBLIVIOUSROUTING_RAECKE_MST_SOLVER_H

#include "raecke_random_mst.h"
#include "raecke_mst_transform.h"
#include "../raecke_framework.h"


/* * This class implements Räcke’s oblivious‐routing solver
 * as in his STOC 2008 paper using the uniform MST algorithm.
 */
// ---------------------- raecke_mst_solver.h ----------------------


#include "raecke_random_mst.h"
#include "raecke_mst_transform.h"
#include "../../solver/solver.h"

class RaeckeMSTSolver : public ObliviousRoutingSolver {
    RaeckeMST mst_algo;
    RaeckeMSTTransform transform;

public:

    RaeckeMSTSolver() = default;
    ~RaeckeMSTSolver() = default;

    void runSolve(const IGraph &g_) override {
        auto g = dynamic_cast<const Graph&>(g_); // cast to Graph
        mst_algo.setGraph(g);
        mst_algo.run(transform);  // pass transform directly
        iteration_count = mst_algo.getIterationCount();
        oracle_running_times = mst_algo.oracle_running_times;
        scaleDownFlow(); // normalize after building
    }

    void storeFlow() override {} // handled during run()

    void scaleDownFlow() {
        std::unordered_map<std::pair<int,int>, double> total_outflow;
        auto routing = transform.getRoutingTable();
        for (const auto& [edge, dmap] : routing) {
            for (const auto& [dem, flow] : dmap) {
                if (edge.first == dem.first ||
                    edge.second == dem.first)
                    total_outflow[dem] += std::abs(flow);
            }
        }

        for (auto& [edge, dmap] : routing) {
            for (auto& [dem, flow] : dmap)
                f_e_st[edge][dem] = flow / total_outflow[dem];
        }
        //f_e_st = routing; // final flow table
    }
};


#endif //OBLIVIOUSROUTING_RAECKE_MST_SOLVER_H