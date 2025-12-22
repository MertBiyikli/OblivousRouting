//
// Created by Mert Biyikli on 24.11.25.
//

#ifndef OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_TRANSFORM_H
#define OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_TRANSFORM_H
#include <vector>
#include <set>
#include <optional>

#include "raecke_oracle_iteration.h"
#include "ckr/ckr_tree_decomposer.h"
#include "../solver/routing_table.h"


/*
 * The efficient Raecke transform class
 * implements the transformation from tree-based routing to
 * an efficient routing table
 */
class IRaeckeTransform {
public:
    // members
    IGraph& g;
    std::vector<OracleTreeIteration>& iterations;
    virtual ~IRaeckeTransform() = default;

    IRaeckeTransform(IGraph& _g, std::vector<OracleTreeIteration>& _iter):g(_g), iterations(_iter) {};

    virtual void distributeDemands(const std::shared_ptr<ITreeNode> &node, double lambda, const std::vector<double>& distance) = 0;

    virtual void removeCycles() = 0;
    virtual std::optional<std::vector<std::pair<int, int>>> findCycleRec(
        int u,
        std::set<int>& analyzed,
        std::vector<int>& stack,
        const std::pair<int, int>& d) = 0;

    virtual std::vector<std::pair<int, int>> findCycle(const std::pair<int, int>& d) = 0;
    void transform();

    virtual bool routingTableIsValid() = 0;
protected:


    virtual void addTree(OracleTreeIteration &iteration);
    std::set<int> collectSubtreeVertices(const std::shared_ptr<ITreeNode> &node);
};



class EfficientRaeckeTransform : public IRaeckeTransform {
protected:
    int root = 0;
    AllPairRoutingTable table;

public:
    virtual ~EfficientRaeckeTransform() = default;

    std::vector<std::unordered_map<std::pair<int, int>, double, PairHash>> edge2demand2flow;
    EfficientRaeckeTransform(IGraph& _g, std::vector<OracleTreeIteration>& _iter):IRaeckeTransform(_g, _iter) {
        table.init((g));
        edge2demand2flow.resize(g.getNumEdges());
    };

    virtual bool routingTableIsValid() override;
    const AllPairRoutingTable& getRoutingTable() const;

    void distributeDemands(const std::shared_ptr<ITreeNode> &node, double lambda, const std::vector<double>& distance) override;

    virtual void removeCycles() override;
    virtual std::optional<std::vector<std::pair<int, int>>> findCycleRec(
        int u,
        std::set<int>& analyzed,
        std::vector<int>& stack,
        const std::pair<int, int>& d) override;

    virtual std::vector<std::pair<int, int>> findCycle(const std::pair<int, int>& d) override;

};



class LinearEfficientRaeckeTransform : public IRaeckeTransform{
protected:
    int root = 0;

public:
    virtual ~LinearEfficientRaeckeTransform() = default;

    LinearEfficientRaeckeTransform(IGraph& _g, std::vector<OracleTreeIteration>& _iter):IRaeckeTransform(_g, _iter){
        linear_tb.init(g);
    };

    LinearRoutingTable linear_tb;

    const LinearRoutingTable& getRoutingTable() const;

    virtual bool routingTableIsValid() override {
        return linear_tb.isValid(g);
    }
private:

    void distributeDemands(const std::shared_ptr<ITreeNode> &node, double lambda, const std::vector<double>& distance) override;

    virtual void removeCycles() override;
    virtual std::optional<std::vector<std::pair<int, int>>> findCycleRec(
        int u,
        std::set<int>& analyzed,
        std::vector<int>& stack,
        const std::pair<int, int>& d) override;

    virtual std::vector<std::pair<int, int>> findCycle(const std::pair<int, int>& d) override;
};


#endif //OBLIVIOUSROUTING_EFFICIENT_RAECKE_CKR_TRANSFORM_H