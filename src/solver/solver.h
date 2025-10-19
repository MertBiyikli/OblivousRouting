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
    bool debug = false;
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


    // TODO: make the adjacency list representation global for all solvers, except maybe the Applegate and Cohen one...
    std::vector<std::vector<int>> adj_f_e_u_id; // per-adjacency list version of f_e_u (stores edge ids)
    std::vector<std::vector<double>> adj_f_e_u; // per-adjacency list version of f_e_u
    void addFlow(int edge_id, int u, double flow) {
        // keep the edges and flows sorted by id
        if (adj_f_e_u_id[edge_id].empty() || adj_f_e_u_id[edge_id].back() < u) {
            adj_f_e_u_id[edge_id].push_back(u);
            adj_f_e_u[edge_id].push_back(flow);
        } else {
            auto it = std::ranges::lower_bound(adj_f_e_u_id[edge_id], u);
            int idx = std::distance(adj_f_e_u_id[edge_id].begin(), it);
            if (it != adj_f_e_u_id[edge_id].end() && *it == u) {
                adj_f_e_u[edge_id][idx] += flow; // accumulate
            } else {
                adj_f_e_u_id[edge_id].insert(it, u);
                adj_f_e_u[edge_id].insert(adj_f_e_u[edge_id].begin() + idx, flow);
            }
        }
    }
};

#endif //OBLIVOUSROUTING_SOLVER_H
