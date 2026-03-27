//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_POINTER_HST_H
#define OBLIVIOUSROUTING_POINTER_HST_H

#include <vector>
#include <memory>

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

bool same_members(const std::vector<int>& a, const std::vector<int>& b);
int calculatePointerHSTHeight(const std::shared_ptr<HSTNode>& node);

#endif //OBLIVIOUSROUTING_POINTER_HST_H