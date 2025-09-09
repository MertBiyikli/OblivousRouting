//
// Created by Mert Biyikli on 05.06.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_TREE_DECOMP_H
#define OBLIVIOUSROUTING_RAECKE_TREE_DECOMP_H


#include <memory>
#include <vector>
#include <unordered_map>
#include "mcct/mcct_derandomized_weighted_solver.h"
#include "../solver/solver.h"



/**
 * This class implements Räcke’s oblivious‐routing solver
 * as in his STOC 2008 paper using the FRT algorithm.
 */

class RaeckeFRT {

    MCCT_Solver m_mcct;
    Graph m_graph;
    double m_lambdaSum;
    std::vector<double> m_lambdas;
    std::vector<FRT_Tree> m_trees;
    std::vector<Graph> m_graphs;
    std::vector< std::unordered_map<std::pair<int, int>, double> > m_idTree2edge2rload;
public:
    bool debug = false;
    // Getters / Setters
    Graph getGraph() const;
    void setGraph(const Graph& g);
    std::vector<FRT_Tree> getTrees() const;
    std::vector<double> getLambdas() const;
    std::vector<Graph> getGraphs() const;

    void run();

    double iterate(int treeIndex);

    FRT_Tree getTree(std::shared_ptr<Graph>& g);
    void computeRLoads(int treeIndex,
                       FRT_Tree& _t,
                       std::shared_ptr<Graph>& copyGraph);
    double getMaxRload(int treeIndex,FRT_Tree& _t);

    void setRequirements(const std::shared_ptr<Graph>& g);
    void computeNewDistances(std::shared_ptr<Graph>& g);
    void normalizeDistance(std::shared_ptr<Graph>& _g, std::unordered_map<std::pair<int, int>, double>& edge2scaledDist);
    double getRloadAllEdges(const std::shared_ptr<Graph>& g) const;

    double getLambdaSum() const {
        return m_lambdaSum;
    }
};





#endif //OBLIVIOUSROUTING_RAECKE_TREE_DECOMP_H
