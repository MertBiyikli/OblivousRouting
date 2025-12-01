//
// Created by Mert Biyikli on 24.11.25.
//

#ifndef OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_TRANSFORM_H
#define OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_TRANSFORM_H
#include <vector>
#include <set>
#include <optional>
#include "../../../datastructures/graph_csr.h"
#include "../../raecke_oracle_iteration.h"
#include "../ckr_tree_decomposer.h"
/*
 * The efficient Raecke CKR Transform class should
 */

class EfficientTreeRoutingTable {
public:
    void init(const int numEdges);
    double getFraction(int e , int s, int t);
    void addFraction(const int e, const int s, const int t, const double fraction);

    std::vector<double> operator[](const int e);

    // store the flows for each commodity
    std::vector<std::vector<std::pair<int, int>>> adj_ids; // adj_ids[e] = [s1, s2, ...] list of commodities for edge e
    std::vector<std::vector<double>> adj_vals; // adj_vals[e] = [f1, f2, ...] list of flows for edge e corresponding to adj_ids
};



class EfficientRaeckeCKRTransform {
protected:
    EfficientTreeRoutingTable table;
    const Graph_csr* g = nullptr;
    std::vector<OracleTreeIteration<std::shared_ptr<TreeNode>, std::vector<double> >>* iterations = nullptr;
public:
    virtual ~EfficientRaeckeCKRTransform() = default;

    void init(IGraph& graph, std::vector<OracleTreeIteration<std::shared_ptr<TreeNode>, std::vector<double> >>& _iters);


    void transform();
    virtual EfficientTreeRoutingTable& addTree(OracleTreeIteration<std::shared_ptr<TreeNode>, std::vector<double> >& iteration);


    const EfficientTreeRoutingTable& getRoutingTable() const;

private:
    void setGraph(Graph_csr& graph);
    void setIterations(std::vector<OracleTreeIteration<std::shared_ptr<TreeNode>, std::vector<double> >>& iterations);
    void distributeDemands(const std::shared_ptr<TreeNode> &node, double lambda, const std::vector<double>& distance);

    std::set<int> collectSubtreeVertices(const std::shared_ptr<TreeNode>& node);
    void normalizeOldSolutionBasedOnNewLambda(double lambda);

    void removeCycles();
    std::optional<std::vector<std::pair<int, int>>> findCycleRec(
        int u,
        std::set<int>& analyzed,
        std::vector<int>& stack,
        const std::pair<int, int>& d);

    std::vector<std::pair<int, int>> findCycle(const std::pair<int, int>& d);

};


#endif //OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_TRANSFORM_H