//
// Created by Mert Biyikli on 24.11.25.
//

#ifndef OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_TRANSFORM_H
#define OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_TRANSFORM_H
#include <vector>
#include <set>
#include <optional>
#include "../../../datastructures/GraphCSR.h"
#include "../../raecke_oracle_iteration.h"
#include "../ckr_tree_decomposer.h"
#include "../../../solver/routing_table.h"


/*
 * The efficient Raecke CKR Transform class should
 */
class EfficientRaeckeCKRTransform {
protected:
    EfficientRoutingTable table;
    const GraphCSR* g = nullptr;
    std::vector<OracleTreeIteration<std::shared_ptr<TreeNode>, std::vector<double> >>* iterations = nullptr;
public:
    virtual ~EfficientRaeckeCKRTransform() = default;

    void init(IGraph& graph, std::vector<OracleTreeIteration<std::shared_ptr<TreeNode>, std::vector<double> >>& _iters);


    void transform();
    virtual EfficientRoutingTable& addTree(OracleTreeIteration<std::shared_ptr<TreeNode>, std::vector<double> >& iteration);


    const EfficientRoutingTable& getRoutingTable() const;

private:
    void setGraph(GraphCSR& graph);
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