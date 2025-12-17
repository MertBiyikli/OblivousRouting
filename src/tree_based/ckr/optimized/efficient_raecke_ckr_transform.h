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
class EfficientRaeckeTransform {
protected:
    int root = 0;
    EfficientRoutingTable table;
    IGraph& g;
    std::vector<OracleTreeIteration>& iterations;
public:
    virtual ~EfficientRaeckeTransform() = default;

    EfficientRaeckeTransform(IGraph& _g, std::vector<OracleTreeIteration>& _iter):g(_g), iterations(_iter) {
        table.init((g.getNumEdges()));
        linear_tb.init(g);
    };

    LinearRoutingTable linear_tb;

    void transform();
    virtual void addTree(OracleTreeIteration &iteration);


    const EfficientRoutingTable& getRoutingTable() const;

private:

    void distributeDemands(const std::shared_ptr<ITreeNode> &node, double lambda, const std::vector<double>& distance);

    std::set<int> collectSubtreeVertices(const std::shared_ptr<ITreeNode> &node);
    void normalizeOldSolutionBasedOnNewLambda(double lambda);

    void removeCycles();
    std::optional<std::vector<std::pair<int, int>>> findCycleRec(
        int u,
        std::set<int>& analyzed,
        std::vector<int>& stack,
        const std::pair<int, int>& d);

    std::vector<std::pair<int, int>> findCycle(const std::pair<int, int>& d);

};



class LinearEfficientRaeckeTransform {
protected:
    int root = 0;
    IGraph& g;
    std::vector<OracleTreeIteration>& iterations;
public:
    virtual ~LinearEfficientRaeckeTransform() = default;

    LinearEfficientRaeckeTransform(IGraph& _g, std::vector<OracleTreeIteration>& _iter):g(_g), iterations(_iter) {
        linear_tb.init(g);
    };

    LinearRoutingTable linear_tb;

    void transform();
    virtual void addTree(OracleTreeIteration &iteration);


    const LinearRoutingTable& getRoutingTable() const;

private:

    void distributeDemands(const std::shared_ptr<ITreeNode> &node, double lambda, const std::vector<double>& distance);

    std::set<int> collectSubtreeVertices(const std::shared_ptr<ITreeNode> &node);
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