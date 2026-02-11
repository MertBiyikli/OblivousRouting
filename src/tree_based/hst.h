//
// Created by Mert Biyikli on 09.02.26.
//

#ifndef OBLIVIOUSROUTING_HST_H
#define OBLIVIOUSROUTING_HST_H

#include <vector>

class HSTNode {
public:
    int center = -1; // center vertex of the cluster
    int id; // unique identifier for the node
    std::vector<std::shared_ptr<HSTNode>> children;
    std::weak_ptr<HSTNode> parent;
    std::vector<int> members; // vertices in the original graph that belong to this node
    explicit HSTNode(int _id) : id(_id) {}
    ~HSTNode() = default;

    const std::vector<std::shared_ptr<HSTNode>>& getChildren() const {return children;}
    const std::vector<int>& getMembers() const {return members;}
};

class HSTLevel {
public:
    double R; // radius of the clusters at this level
    std::vector<int> centers;
    std::vector<int> owner;
};
#endif //OBLIVIOUSROUTING_HST_H