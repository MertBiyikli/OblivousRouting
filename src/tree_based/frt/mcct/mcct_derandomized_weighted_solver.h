//
// Created by Mert Biyikli on 20.05.25.
//

#ifndef OBLIVIOUSROUTING_MCCT_DERANDOMIZED_WEIGHTED_SOLVER_H
#define OBLIVIOUSROUTING_MCCT_DERANDOMIZED_WEIGHTED_SOLVER_H


#include "../../../datastructures/graph.h"
#include "Tree.h"

#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include <memory>
#include <cmath>
#include <iostream>
#include <unordered_map>


class MCCT_Solver {
    bool debug = false;
    FRT_Tree tree;
    std::shared_ptr<Graph> graph;
    std::unordered_set<double> betas;
    double bestBeta;
    std::vector<int> verticesPermutation;
    std::unordered_map<int, std::unordered_map<int, double>> idVertex2idVertex2demand; // Set of demands
    std::map<std::pair<int, int>, int> demand2levelIncluded; // Map of demands to their levels in the tree
    public:


    FRT_Tree getBestTree(bool debug = false);
    void computeBestTree(bool debug = false);

    void computeBestBetaAndPermutation();
    void computeBetas(bool debug = false);
    double computeExpectation(double beta, std::vector<int>& verticesPermutation, std::set<int>& unsettledVertices, bool debug = false);


    void addDemand(int v1, int v2, double demand) {
        if(v1 < 0 || v2 < 0) {
            throw std::invalid_argument("Vertex IDs must be non-negative.");
        }
        if (v1 == v2) {
            throw std::invalid_argument("Cannot add demand from a vertex to itself.");
        }
        // Ensure that v1 and v2 are in the idVertex2idVertex2demand map
        auto& it = idVertex2idVertex2demand[v1];
        if (it.empty()) {
            idVertex2idVertex2demand[v1] = std::unordered_map<int, double>();
        }

        it[v2] = demand;
    }

    void removeDemands(double beta, int v);

    FRT_Tree getTree() const {
        return tree;
    }

    void setTree(const FRT_Tree& t) {
        tree = t;
    }

    Graph& getGraph() const {
        if(!graph) {
            throw std::runtime_error("Graph is not set.");
        }
        return *graph;
    }

    void setGraph(Graph& g) {
        graph = std::make_shared<Graph>(g);
    }

    void reset() {
        betas.clear();
        bestBeta = 0.0;
        verticesPermutation.clear();
        idVertex2idVertex2demand.clear();
        demand2levelIncluded.clear();
    }

    const std::unordered_set<double>& getBetas() const {
        return betas;
    }

    void setBetas(const std::unordered_set<double>& b) {
        betas = b;
    }

    double getBestBeta() const {
        return bestBeta;
    }

    void setBestBeta(double b) {
        bestBeta = b;
    }

    const std::vector<int>& getVerticesPermutation() const {
        return verticesPermutation;
    }

    void setVerticesPermutation(const std::vector<int>& perm) {
        verticesPermutation = perm;
    }
};



/*
class MCCTDerandomizedWeightedSolver {
public: // TODO: change this to private
    Tree tree;
    std::shared_ptr<Graph> graph;
    std::set<double> betas;
    double bestB;
    std::vector<int> verticesPermutation;
    std::map<int, std::map<int, double>> idVertex2idVertex2demand; // Set of demands
    std::map<std::pair<int, int>, int> demand2levelIncluded;

    // After you build the tree, for each tree‐edge (u→v) you will
    // insert {u,v}→arcPtr into this map.
    std::map<std::pair<int,int>, std::shared_ptr<Edge> > treeArcMap_;

public:
    MCCTDerandomizedWeightedSolver();

    // Return the std::shared_ptr<Edge>  that MCCT used to connect centerU → centerV,
    // or nullptr if it doesn’t exist.
    std::shared_ptr<Edge>  getArcBetween(int centerU, int centerV) const {
        auto it = treeArcMap_.find({centerU, centerV});
        return (it == treeArcMap_.end() ? nullptr : it->second);
    }

    // Add demand from vertex v1 to vertex v2
    void addDemand(int v1, int v2, double demand);

    // Get the tree
    Tree getTree() const;

    // Set the tree
    void setTree(const Tree& tree);

    // Get the graph
    std::shared_ptr<Graph> getGraph() const;

    // Set the graph
    void setGraph(std::shared_ptr<Graph> graph);

    // Get betas
    std::set<double> getBetas() const;

    // Set betas
    void setBetas(const std::set<double>& betas);

    // Get the best beta value
    double getBestB() const;

    // Set the best beta value
    void setBestB(double bestB);

    // Get the vertices permutation
    std::vector<int> getVerticesPermutation() const;

    // Set the vertices permutation
    void setVerticesPermutation(const std::vector<int>& verticesPermutation);

    // Find best beta and vertex permutation
    void findBestBetaAndPermutation();

    // Remove demands whose cost is already settled
    void removeDemands(double beta, int bestV);

    // Compute the expected cost of the decomposition tree given a certain beta
    double computeExpectation(double beta, std::vector<int>& verticesPermutation, std::set<int>& unsettledVertices, bool debug = false);

    // Compute all meaningful betas
    void computeBetas();

    // Get the best tree after computation
    Tree getBestTree();

    // Compute the best tree using the best beta and permutation
    void computeBestTree();

    // Reset solver state
    void reset();
};
*/


#endif //OBLIVIOUSROUTING_MCCT_DERANDOMIZED_WEIGHTED_SOLVER_H
