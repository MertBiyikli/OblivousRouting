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
#include "../graph.h"
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
        double totalRLoadsAllEdges = this->getRloadAllEdges(g);
        std::unordered_map<std::pair<int, int>, double> edge2scaledDist;

        for(int u = 0; u<g.getNumNodes(); u++) {
            for (const auto& v : g.neighbors(u)) {
                double totalrLoad = 0.0;
                for (size_t i = 0; i < m_trees.size(); ++i) {
                    auto it = this->m_idTree2edge2rload[i].find({u, v});
                    double rLoad = (it != this->m_idTree2edge2rload[i].end() ? it->second : 0.0);
                    totalrLoad += rLoad * m_lambdas[i];
                }
                double num         = std::exp(totalrLoad) / g.getEdgeCapacity(u, v);
                double newDistance = num / totalRLoadsAllEdges;

                if(std::isinf(newDistance) || std::isinf(num)) {
                    throw std::runtime_error("Infinity encountered in newDistance or Exponential calculation for edge: " + std::to_string(u) + " → " + std::to_string(v));
                }

                if(std::isnan(newDistance)) {
                    throw std::runtime_error("NaN encountered in newDistance calculation for edge: " + std::to_string(u) + " → " + std::to_string(v));
                }

                edge2scaledDist[{u, v}] = newDistance;
            }
        }

        // Normalize distances
        normalizeDistance(g, edge2scaledDist);
    }


    virtual void normalizeDistance(Graph &_g, std::unordered_map<std::pair<int, int>, double> &edge2scaledDist) {
        double minDistance = std::numeric_limits<double>::infinity();
        for(int u = 0; u < _g.getNumNodes(); u++) {
            for (const auto& v : _g.neighbors(u)) {
                double scaledDist = edge2scaledDist[{u, v}];
                if (scaledDist < minDistance) {
                    minDistance = scaledDist;
                }
            }
        }

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