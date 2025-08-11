//
// Created by Mert Biyikli on 05.06.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_TREE_DECOMP_H
#define OBLIVIOUSROUTING_RAECKE_TREE_DECOMP_H


#include <memory>
#include <vector>
#include <unordered_map>
#include <map>
#include <queue>
#include <algorithm>
#include <cmath>
#include "../utils/hash.h"
#include "mcct/mcct_derandomized_weighted_solver.h"
#include "../solver/solver.h"

class RaeckeFRT {
    bool debug = false;
    MCCT_Solver m_mcct;
    Graph m_graph;
    double m_lambdaSum;
    std::vector<double> m_lambdas;
    std::vector<FRT_Tree> m_trees;
    std::vector<Graph> m_graphs;
    std::vector< std::unordered_map<std::pair<int, int>, double> > m_idTree2edge2rload;
public:

    // Getters / Setters
    Graph getGraph() const;
    void setGraph(const Graph& g);

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


    std::vector<FRT_Tree> getTrees() const {
        return m_trees;
    }

    std::vector<double> getLambdas() const {
        return m_lambdas;
    }

    std::vector<Graph> getGraphs() const {
        return m_graphs;
    }

    double getLambdaSum() const {
        return m_lambdaSum;
    }
};

using EdgeWeights = std::unordered_map<std::shared_ptr<Edge>, double>;

/**
 * This class implements Räcke’s oblivious‐routing solver
 * as in his STOC 2008 paper using the FRT algorithm.
 */

/*
class RackeObliviousRoutingSolver : public OblviviousRoutingSolver {
public:
    RackeObliviousRoutingSolver()
        : graph_(nullptr){
        mcct_ = std::make_unique<MCCTDerandomizedWeightedSolver>();
    };

    virtual ~RackeObliviousRoutingSolver();

    virtual void solve(const Graph& graph) override;

    // The main entry‐point
    void createTreesAndLambda();

    // Getters / Setters
    std::vector<std::shared_ptr<Graph>>&        getGraphs();
    void                                        setGraph(const std::shared_ptr<Graph>& g);
    std::vector<double>                         getLambdas() const;
    std::vector<std::shared_ptr<Tree>>          getTrees() const;
    std::shared_ptr<Tree>                       getTree(std::shared_ptr<Graph>& g);



    void printObliviousRoutingTable() const;


private:
    // Helper methods to match the Java version
    void setRequirements(const std::shared_ptr<Graph>& g);
    void computeRLoads(int treeIndex,
                       const std::shared_ptr<Tree>& t,
                       const std::shared_ptr<Graph>& copyGraph);


    double getMaxRLoad(int treeIndex,
                           const std::shared_ptr<Tree>& t);


    EdgeWeights computeNewDistances(std::shared_ptr<Graph>& g);


    void normalizeDistance(std::shared_ptr<Graph>& _g, EdgeWeights& arc2scaledDist);

    double getRloadAllEdges(const std::shared_ptr<Graph>& g) const;
    std::vector<EdgeWeights> GetEdge2Load() const;


    std::shared_ptr<Graph>                               graph_;
    std::vector<std::shared_ptr<Tree>>                    trees_;
    std::vector<double>                                   lambdas_;
    std::vector<std::shared_ptr<Graph>>                   graphs_;
    // For each tree‐index, and for each Arc*, store its “R‐load”:
    std::vector<EdgeWeights>         idTree2arc2rLoad_; // todo: change to Arc Ids instead of pointers
    std::unique_ptr<MCCTDerandomizedWeightedSolver>      mcct_;
};
*/



#endif //OBLIVIOUSROUTING_RAECKE_TREE_DECOMP_H
