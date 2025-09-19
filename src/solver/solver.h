//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_SOLVER_H
#define OBLIVOUSROUTING_SOLVER_H

#include <vector>
#include "../graph.h"
#include "../utils/hash.h"

class ObliviousRoutingSolver { // TODO: rename the class name
public:
    std::vector<double> oracle_running_times;
    int iteration_count = 0;
    std::unordered_map<std::pair<int, int> , std::unordered_map<std::pair<int, int>, double >> f_e_st;
    std::vector<std::vector<double>> m_routingTable;
    ObliviousRoutingSolver() = default;
    virtual ~ObliviousRoutingSolver() = default;

    // ToDo: Think of an efficient way of storing the routing table
    virtual void solve(const Graph& graph) = 0;
    virtual void storeFlow() = 0;


    int GetIterationCount() const {
        return iteration_count;
    }

    void printFlow() const {
        for (const auto& [edge, flow_map] : f_e_st) {
            std::cout << "Edge (" << edge.first << ", " << edge.second << "): ";
            for (const auto& [commodity, value] : flow_map) {
                std::cout << "  Commodity (" << commodity.first << ", " << commodity.second << ") -> Flow: " << value << "; ";
            }
            std::cout << std::endl;
        }
    }
};

#endif //OBLIVOUSROUTING_SOLVER_H
