//
// Created by Mert Biyikli on 19.09.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_BASE_H
#define OBLIVIOUSROUTING_RAECKE_BASE_H
#include <chrono>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <stdexcept>
#include "../datastructures/graph.h"
#include "../utils/hash.h"


/** This is the base class for Raecke's framework
 * It is templated by the type of tree used in the decomposition
 * (e.g., FRT_Tree or MSTTree)
 */
template<typename Tree>
class RaeckeBase {
protected:
    double m_lambdaSum;
    Graph m_graph;

    std::vector<double> m_lambdas;
    std::vector<Graph> m_graphs;

    std::vector<Tree> m_trees;
    std::vector< std::unordered_map<std::pair<int, int>, double> > m_idTree2edge2rload;

public:
    std::vector<double> oracle_running_times; // TODO: remove this
    std::vector<double> pure_oracle_running_times; // TODO: remove this
    virtual void init(Graph& g) = 0;
    virtual Tree getTree(Graph& g) = 0;
    virtual void computeRLoads(int treeIndex,
                       Tree& _t,
                       Graph& copyGraph) = 0;


    void run() {
        init(m_graph);
        int id = 0;
        m_lambdaSum = 0.0;
        while (m_lambdaSum < 1.0) {
            auto start = std::chrono::high_resolution_clock::now();
            m_lambdaSum += iterate(id);

            oracle_running_times.push_back((std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start)).count());
            id++;
        }
    }

    virtual double iterate(int treeIndex) {
        // m_graph.print();
        Tree t = getTree(m_graph);
        m_trees.push_back(t);

        computeRLoads(treeIndex, t, m_graph);
        double l = getMaxRload(treeIndex);
        double delta = std::min(1/l, 1-m_lambdaSum);
        m_lambdas.push_back(delta);

        m_graphs.push_back(Graph(m_graph));
        return delta;
    }




    Graph getGraph() const {return m_graph;}
    void setGraph(const Graph& g) { m_graph = g; }
    const std::vector<double>& getLambdas() const { return m_lambdas; }
    const std::vector<Graph>& getGraphs() const { return m_graphs; }
    std::vector<Graph> getGraphs() { return m_graphs; }
    const std::vector<double>& getOracleTimes() const { return oracle_running_times; }
    const std::vector<Tree>& getTrees() const { return m_trees; }
    std::vector<Tree> getTrees() { return m_trees; }
    double getLambdaSum() const { return m_lambdaSum; }



    virtual double getMaxRload(int treeIndex) {
        double maxRatio = 0.0;
        for(const auto& [edge, rLoad] : m_idTree2edge2rload[treeIndex]) {
            if (rLoad > maxRatio) {
                maxRatio = rLoad;
            }
        }
        return maxRatio;
    }


    virtual double getRloadAllEdges(const Graph& g) const{
        double totalRLoadAllEdges = 0.0;
        for(int u = 0; u < g.getNumNodes(); u++) {
            for (const auto& v : g.neighbors(u)) {
                double totalPerEdge = 0.0;

                std::pair<int, int> edge = {u, v};
                // Sum over all trees i:  rLoad(i,arc) * lambda[i]
                for (size_t i = 0; i < m_trees.size(); ++i) {
                    auto it = this->m_idTree2edge2rload[i].find(edge);
                    if(it != this->m_idTree2edge2rload[i].end()) {
                        totalPerEdge += it->second * m_lambdas[i];
                    } else {
                        // If the edge is not found, we assume the r-load is 0
                        totalPerEdge += 0.0;
                    }
                }

                // Add e^(totalPerEdge) to the grand total
                totalRLoadAllEdges += std::exp(totalPerEdge);
            }
        }
        return totalRLoadAllEdges;
    }


    virtual void computeNewDistances(Graph &g) {
        constexpr double EPS = 1e-12;

        // --- 1) Compute totalrLoad for every edge and track max (for log-sum-exp) ---
        std::unordered_map<std::pair<int,int>, double> total_r_by_edge;
        total_r_by_edge.reserve(static_cast<size_t>(g.getNumEdges()) * 2);

        double max_total_r = -std::numeric_limits<double>::infinity();

        for (int u = 0; u < g.getNumNodes(); ++u) {
            for (int v : g.neighbors(u)) {
                double totalrLoad = 0.0;
                for (size_t i = 0; i < m_trees.size(); ++i) {
                    auto it = m_idTree2edge2rload[i].find({u, v});
                    if (it != m_idTree2edge2rload[i].end()) {
                        totalrLoad += it->second * m_lambdas[i];
                    }
                }
                total_r_by_edge[{u, v}] = totalrLoad;
                if (totalrLoad > max_total_r) max_total_r = totalrLoad;
            }
        }

        // If no trees yet (first iteration), we still need (finite) distances.
        if (!std::isfinite(max_total_r)) max_total_r = 0.0;

        // --- 2) Compute denominator using the same shift (log-sum-exp) ---
        // denom = Σ_e exp(totalrLoad(e))
        //       = exp(M) * Σ_e exp(totalrLoad(e) - M)
        // We only need the scaled sum: sumExp = Σ_e exp(totalrLoad(e) - M)
        double sumExp = 0.0;
        for (int u = 0; u < g.getNumNodes(); ++u) {
            for (int v : g.neighbors(u)) {
                double tr = total_r_by_edge[{u, v}] - max_total_r;
                // exp of (possibly negative) number; safe
                sumExp += std::exp(tr);
            }
        }
        // Guard against degenerate graphs
        if (sumExp <= 0.0 || !std::isfinite(sumExp)) {
            throw std::runtime_error("computeNewDistances: invalid sumExp in log-sum-exp normalization.");
        }

        // --- 3) Update distances with the same shift in numerator & denominator ---
        // newDistance(u,v) = [exp(totalrLoad(u,v)) / cap(u,v)] / Σ_e exp(totalrLoad(e))
        //                  = [exp(totalrLoad(u,v) - M) / cap(u,v)] / Σ_e exp(totalrLoad(e) - M)
        //                  = safe, because both numerator and denominator are well-scaled.
        std::unordered_map<std::pair<int,int>, double> edge2scaledDist;
        edge2scaledDist.reserve(total_r_by_edge.size());

        for (int u = 0; u < g.getNumNodes(); ++u) {
            for (int v : g.neighbors(u)) {
                double cap = g.getEdgeCapacity(u, v);
                if (cap < EPS) cap = EPS;

                double tr_shifted = total_r_by_edge[{u, v}] - max_total_r;
                double num_scaled = std::exp(tr_shifted) / cap;  // safe

                if (!std::isfinite(num_scaled)) {
                    throw std::runtime_error("computeNewDistances: non-finite num_scaled for edge (" +
                                             std::to_string(u) + "," + std::to_string(v) + ").");
                }

                double newDistance = num_scaled / sumExp;        // also safe
                if (!std::isfinite(newDistance)) {
                    throw std::runtime_error("computeNewDistances: non-finite newDistance for edge (" +
                                             std::to_string(u) + "," + std::to_string(v) + ").");
                }

                edge2scaledDist[{u, v}] = newDistance;
            }
        }

        // --- 4) Normalize distances to keep min distance = 1 (your original policy) ---
        normalizeDistance(g, edge2scaledDist);
    }


    virtual void normalizeDistance(Graph &_g, std::unordered_map<std::pair<int, int>, double> &edge2scaledDist) {
        double minDistance = std::numeric_limits<double>::infinity();
        for(int u = 0; u < _g.getNumNodes(); u++) {
            for (const auto& v : _g.neighbors(u)) {
                auto it = edge2scaledDist.find({u, v});
                if (it == edge2scaledDist.end()) continue;
                double d = it->second;
                if (std::isfinite(d) && d < minDistance) minDistance = d;
            }
        }
        if (!std::isfinite(minDistance) || minDistance <= 0.0) minDistance = 1.0;
        for(int u = 0; u < _g.getNumNodes(); u++) {
            for (const auto &v: _g.neighbors(u)) {
                double arc2scaledDistValue = edge2scaledDist[{u, v}];
                double newDistance = arc2scaledDistValue / minDistance;

                if (std::isnan(newDistance)) {
                    throw std::runtime_error(
                            "NaN encountered in newDistance for edge: " + std::to_string(u) + " → " + std::to_string(v));
                }
                if (newDistance < 1) {
                    newDistance = 1.0; // Ensure minimum distance is 1}
                }
                _g.updateEdgeDistance(u, v, newDistance);
            }
        }
    }
};

#endif //OBLIVIOUSROUTING_RAECKE_BASE_H