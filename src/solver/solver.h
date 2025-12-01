//
// Created by Mert Biyikli on 11.05.25.
//

#ifndef OBLIVOUSROUTING_SOLVER_H
#define OBLIVOUSROUTING_SOLVER_H

#include <vector>
#include "../datastructures/graph.h"
#include "../utils/hash.h"
#include <unordered_map>
//#include <boost/unordered_map.hpp>

class Graph_csr;



class ObliviousRoutingSolver { // TODO: rename the class name
protected:
    std::unique_ptr<IGraph> graph;
public:
    bool debug = false;
    std::vector<double> oracle_running_times;
    std::vector<double> pure_oracle_running_times;
    int iteration_count = 0;
    std::unordered_map<std::pair<int, int> , std::unordered_map<std::pair<int, int>, double, PairHash >, PairHash> f_e_st;

    std::vector<std::vector<double>> m_routingTable;
    ObliviousRoutingSolver() = default;
    virtual ~ObliviousRoutingSolver() = default;


    void solve() {
        if (graph) {
            runSolve(*graph);
            graph->resetEdgeWeights();
        }else {
            throw std::runtime_error("Graph not set in ObliviousRoutingSolver.");
        }
    }

    template<class T>
    T& graphAs() {
        return dynamic_cast<T&>(*graph);
    }

    // ToDo: Think of an efficient way of storing the routing table

    // TODO: use the graph interface to support both CSR and Adjacency list
    virtual void runSolve(const IGraph& graph) = 0;


    template<class GraphType>
    void setGraph(GraphType g) {
        graph = std::make_unique<GraphType>(std::move(g));
    }
    virtual void storeFlow() = 0;


    int GetIterationCount() const {
        return iteration_count;
    }

    void printFlow_() const {
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
    // -----------------------------------------------------------------------------
    // ⚡ Optimized addFlow (manual binary search, cache-friendly)
    // -----------------------------------------------------------------------------
#if defined(__GNUC__) || defined(__clang__)
    __attribute__((always_inline))
    #elif defined(_MSC_VER)
    __forceinline
    #endif
    inline void addFlow(int edge_id, int u, double flow) noexcept {
        auto &ids  = adj_f_e_u_id[edge_id];
        auto &vals = adj_f_e_u[edge_id];
        const size_t len = ids.size();

        // Reserve to avoid reallocations
        if (len == 0) [[unlikely]] {
            ids.reserve(8);
            vals.reserve(8);
            ids.push_back(u);
            vals.push_back(flow);
            return;
        }

        // ✅ Linear scan for small adjacency lists (<= 16 entries)
        if (len <= 16) [[likely]] {
            for (size_t i = 0; i < len; ++i) {
                if (ids[i] == u) {
                    vals[i] += flow;
                    return;
                }
            }
            ids.push_back(u);
            vals.push_back(flow);
            return;
        }

        // ✅ Hand-rolled binary search (faster than std::lower_bound)
        size_t lo = 0, hi = len;
        while (lo < hi) {
            const size_t mid = (lo + hi) >> 1;
            const int mid_val = ids[mid];
            if (mid_val < u)
                lo = mid + 1;
            else
                hi = mid;
        }

        if (lo < len && ids[lo] == u) {
            vals[lo] += flow;
        } else {
            // Usually append to the end — amortized O(1)
            if (lo == len) {
                ids.push_back(u);
                vals.push_back(flow);
            } else {
                ids.insert(ids.begin() + static_cast<long>(lo), u);
                vals.insert(vals.begin() + static_cast<long>(lo), flow);
            }
        }
    }


    void printFlow() {
        for (int e = 0; e < adj_f_e_u.size(); e++) {
            int u = e; // edge index
            for (int idx = 0; idx < adj_f_e_u[e].size(); ++idx) {
                int src = adj_f_e_u_id[e][idx];
                double flow = adj_f_e_u[e][idx];
                std::cout << "Edge " << e << " Commodity (" << src << " -> x_fixed): Flow = " << flow << "\n";
            }
        }
    }

    double getFlowForEdgeCommodity(int edge_id, int s, int t) const {
        for (int idx = 0; idx < adj_f_e_u_id[edge_id].size(); ++idx) {
            if (adj_f_e_u_id[edge_id][idx] == s) {
                return adj_f_e_u[edge_id][idx];
            }
        }
        return 0.0;

    }
};

#endif //OBLIVOUSROUTING_SOLVER_H
