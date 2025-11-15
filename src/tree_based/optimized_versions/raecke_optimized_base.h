//
// Created by Mert Biyikli on 05.11.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_OPTIMIZED_BASE_H
#define OBLIVIOUSROUTING_RAECKE_OPTIMIZED_BASE_H

#include <chrono>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <absl/strings/str_format.h>
#include <absl/strings/internal/str_format/extension.h>

#include "../../graph.h"
#include "../../graph_csr.h"
#include "../../solver/solver.h"
#include "frt/raecke_frt_transform_opt.h"


template<typename Tree, typename Transform>
class RaeckeOptimizedBase : public ObliviousRoutingSolver {
protected:
    double m_lambdaSum;
    Graph_csr m_graph;

    std::vector<double> m_lambdas;
    std::vector<Graph_csr> m_graphs;
    std::vector<Tree> m_trees;
    Transform m_transform;
    FlowStore flow_store;

    // for now keep an extra vector that is managed by the RaeckeOptimizedBase class
    // that is keep the edges of the graphs as pairs sorted in a vector
    // and keep also the capacities, brother...
    // TODO: move this later on in the graph class and embrace CSR format for the base graph type
    /*
    std::vector<std::pair<int, int>> m_edges;
    std::vector<int> head; // store the first appearance of an edge endpoint for faster access
    std::vector<double> m_edge_capacities;
    std::vector<double> m_edge_distances;*/

    std::vector<std::vector<double>> m_tree2edge2rloads; // edge indexed r-loads for each tree

public:


    virtual Tree getTree(Graph_csr& g) = 0;
    virtual void computeRLoads(int treeIndex,
                       Tree& _t,
                       Graph_csr& copyGraph) = 0;


    void solve(const Graph & g) override{
        this->init(g);
        this->run();
        // scale the flow to meet unit flow

        iteration_count = m_lambdas.size();
    }

    void storeFlow() override {

        // Combine all flows accumulated in flow_store into f_e_st
        f_e_st.clear();

        for (int e = 0; e < m_graph.getNumEdges(); ++e) {
            const auto [u, v] = m_graph.edges[e];

            for (const auto& [packed, val] : flow_store.flows_on_edge(e)) {
                auto [s, t] = FlowStore::unpack(packed);
                if (std::abs(val) < 1e-15) continue;

                f_e_st[{u, v}][{s, t}] += val;
            }
        }
        scaleDownFlow();
    }

    void scaleDownFlow() {
        // scale the flow to meet unit flow
        std::unordered_map<std::pair<int, int>, double > outgoingflow_per_commodity;

        for ( const auto& [edge, flowMap]:f_e_st) {
            for (const auto& [com, flow_value]  : flowMap) {
                if ( flow_value < 1e-15 ) continue; // ignore zero flows
                if (!outgoingflow_per_commodity.contains(com) )
                    outgoingflow_per_commodity[com] = 0;

                if (edge.first == com.first
                    || edge.second == com.first) {
                    outgoingflow_per_commodity[com] += std::abs(flow_value);
                    }
            }
        }

        // scale the flow values to meet one unit of flow per commodity
        for ( auto& [edge, flowMap]:f_e_st) {
            for (auto& [com, flow_value] : flowMap) {
                flow_value /= outgoingflow_per_commodity[com];
            }
        }
    }

    void init(const Graph& g) {
        // initialize the base graph
        // convert to csr
        Graph_csr csr_g;
        for (int i = 0; i < g.getNumNodes(); i++) {
            for (const auto& j : g.neighbors(i)) {
                csr_g.addEdge(i, j, g.getEdgeCapacity(i, j));
            }
        }
        csr_g.n = g.getNumNodes(),csr_g.m = g.getNumEdges();
        csr_g.preprocess();
        m_graph = csr_g;
    }


    // base functions
    void run() {
        flow_store.begin(m_graph.getNumEdges());
        int id = 0;
        m_lambdaSum = 0.0;
        while (m_lambdaSum < 1.0) {
            auto start = std::chrono::high_resolution_clock::now();
            m_lambdaSum += iterate(id);

            // right after here already push the computed flow

            oracle_running_times.push_back((std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start)).count());
            id++;
        }
        flow_store.finalize();
    }

    virtual double iterate(int treeIndex) {
        m_graph.print();
        Tree t = getTree(m_graph);
        m_trees.push_back(t);

        computeRLoads(treeIndex, t, m_graph);
        double l = getMaxRload(treeIndex);
        double delta = std::min(1/l, 1-m_lambdaSum);
        m_lambdas.push_back(delta);
        m_graphs.push_back(Graph_csr(m_graph));

        if ( debug) {
            std::cout << m_lambdas.size() << std::endl;
            std::cout << "added tree with max r-load: " << absl::StrFormat("%.6f", l) << ", lambda: " << absl::StrFormat("%.6f", delta) << ", lambda sum: " << absl::StrFormat("%.6f", m_lambdaSum) << std::endl;
        }
        return delta;
    }


    std::vector<std::tuple<int,int,double>> getAllDemands(int max_pairs = -1) const {
        std::vector<std::tuple<int,int,double>> demands;
        const int n = m_graph.getNumNodes();

        for (int s = 0; s < n; ++s) {
            for (int t = s + 1; t < n; ++t) {
                demands.emplace_back(s, t, 1.0);  // unit demand
                if (max_pairs > 0 && (int)demands.size() >= max_pairs)
                    return demands;
            }
        }
        return demands;
    }

    virtual double getMaxRload(int treeIndex) {
        double maxRatio = 0.0;
        for(const auto& rLoad : m_tree2edge2rloads[treeIndex]) {
            if (rLoad > maxRatio) {
                maxRatio = rLoad;
            }
        }
        return maxRatio;
    }
    virtual double getRloadAllEdges() const{
        constexpr double EPS = 1e-12;

        double totalRLoadAllEdges = 0.0;

        for (int e = 0; e < m_graph.edges.size(); e++) {
            double totalRloadEdge = 0.0;
            for (int i_tree = 0; i_tree < m_tree2edge2rloads.size(); i_tree++) {
                const auto& tree_rloads = m_tree2edge2rloads[i_tree];
                if (tree_rloads[e]>EPS) {
                    totalRloadEdge += tree_rloads[e]*m_lambdas[i_tree];
                }
            }
            totalRLoadAllEdges += std::exp(totalRloadEdge);
        }

        return totalRLoadAllEdges;
    }
    virtual void computeNewDistance(Graph_csr& graph) {
        constexpr double EPS = 1e-12;

        const int m = graph.getNumEdges();
        // compute total rload per edge
        std::vector<double> total_rload_per_edge;
        total_rload_per_edge.resize(graph.getNumEdges(), 0.0);
        for(int edge_id = 0 ; edge_id < graph.getNumEdges() ; edge_id++) {
            for (int idx_tree = 0; idx_tree < m_tree2edge2rloads.size(); idx_tree++) {
                const auto& tree_rloads = m_tree2edge2rloads[idx_tree];
                if (tree_rloads[edge_id] < EPS) continue;
                total_rload_per_edge[edge_id] += tree_rloads[edge_id]*m_lambdas[idx_tree];
            }
        }
        // compute max total rload
        double max_total_rload = -std::numeric_limits<double>::infinity();
        for (const auto& rload : total_rload_per_edge) {
            if (rload > max_total_rload) {
                max_total_rload = rload;
            }
        }


        // If no trees yet (first iteration), we still need (finite) distances.
        if (!std::isfinite(max_total_rload)) max_total_rload = 0.0;

        // compute sum exp
        double sumExp = 0.0;
        for (const auto& rload : total_rload_per_edge) {
            double tr = rload - max_total_rload;
            sumExp += std::exp(tr);
        }

        // Guard against degenerate graphs
        if (sumExp <= 0.0 || !std::isfinite(sumExp)) {
            throw std::runtime_error("computeNewDistances: invalid sumExp in log-sum-exp normalization.");
        }

        // --- 3) Compute new distance exactly like the naive version ---
        std::vector<double> new_distances(m);
        for (int e = 0; e < m; ++e) {
            double cap = graph.getEdgeCapacity(e);
            if (cap < EPS) cap = EPS;

            double shifted = total_rload_per_edge[e] - max_total_rload;
            double num = std::exp(shifted) / cap;
            double newDist = num / sumExp;

            if (!std::isfinite(newDist)) {
                auto [u, v] = graph.edges[e];
                throw std::runtime_error(
                    "non-finite distance at edge (" +
                    std::to_string(u) + "," + std::to_string(v) + ")");
            }
            new_distances[e] = newDist;
        }


        // --- 4) normalize the distances
        double min_dist = std::numeric_limits<double>::infinity();
        for (const auto& d : new_distances) {
            // exclude zero distances
            if (d < EPS) continue;
            if (d < min_dist
                && std::isfinite(d)) {
                min_dist = d;
            }
        }
        if (!std::isfinite(min_dist) || min_dist <= 0.0) min_dist = 1.0;

        for (int e = 0; e < m; ++e) {
            double normalized_distance = new_distances[e] / min_dist;
            if (std::isnan(normalized_distance)) {
                throw std::runtime_error(
                        "NaN encountered in newDistance for edge: " + std::to_string(graph.edges[m].first) + " → " + std::to_string(graph.edges[m].second));
            }
            if (normalized_distance < 1) {
                normalized_distance = 1.0; // Ensure minimum distance is 1}
            }
            new_distances[e] = normalized_distance;
        }

        for (int e = 0; e < m; ++e) {
            graph.updateDistance(e, new_distances[e]);
        }

    }
};


#endif //OBLIVIOUSROUTING_RAECKE_OPTIMIZED_BASE_H