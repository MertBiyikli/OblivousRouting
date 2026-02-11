//
// Created by Mert Biyikli on 11.12.25.
//

#ifndef OBLIVIOUSROUTING_EFFICIENT_FRT_H
#define OBLIVIOUSROUTING_EFFICIENT_FRT_H
#include <memory>
#include <unordered_set>
#include "../../datastructures/IGraph.h"
#include "frt_node.h"

class EfficientFRT {
public:
    IGraph& graph;
    std::vector<std::vector<std::shared_ptr<EfficientFRTTreeNode>>>  level2nodes;
    EfficientFRT(IGraph& g) : graph(g), bestBeta(0),vertex2vertex(g.getNumNodes()), vertex2vertex_value(g.getNumNodes())  {
    };
    ~EfficientFRT() = default;

    std::unordered_set<double> betas;
    double bestBeta;
    std::vector<int> verticesPermutation;

    // TODO: this can be optimized !!!!
    std::unordered_map<int, std::unordered_map<int, double>> idVertex2idVertex2demand; // Set of demands

    // instead of nested maps , try to keep an adjacency list representation for the demands
    std::vector<std::vector<int>> vertex2vertex; // demand from u to v: vertex2vertex2demand[u][v]
    std::vector<std::vector<double>> vertex2vertex_value;
    std::unordered_map<std::pair<int, int>, int, PairHash> demand2levelIncluded; // Map of demands to their levels in the tree

    std::shared_ptr<EfficientFRTTreeNode> getTree();

    std::shared_ptr<EfficientFRTTreeNode> getTreeBottomUp();
    void computeCurrentPartition(std::vector<int>& centers, std::vector<int>& clusters, const double scale);

    std::shared_ptr<EfficientFRTTreeNode> getBestTree();
    void updateEdgeDistances(const std::vector<double>& distances);

    void computeBetaAndRandomPermutation();

    void computeBestBetaAndPermutation();
    void computeBetas();

    double computeExpectation(double beta, std::unordered_set<int>& allVertices, const std::vector<int>& currentPermutation);
    void removeDemands(double beta, int v);
    void addDemands(int u, int v, double demand);

    void cleanUpTree(std::shared_ptr<EfficientFRTTreeNode>& node);


};

#endif //OBLIVIOUSROUTING_EFFICIENT_FRT_H