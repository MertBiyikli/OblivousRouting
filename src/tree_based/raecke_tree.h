//
// Created by Mert Biyikli on 09.12.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_TREE_H
#define OBLIVIOUSROUTING_RAECKE_TREE_H
#include <memory>
// Tree interface, abstract base class
class ITreeNode {
public:
    int center = -1;
    virtual ~ITreeNode() = default;
    virtual std::vector<std::shared_ptr<ITreeNode>> getChildren() = 0;
    virtual std::shared_ptr<ITreeNode> getParent() const = 0;
    virtual const std::vector<int>& getMembers() const = 0;

    bool isHierarchical() {
        std::vector<char> visited(getMembers().size(), 0);
        for (const auto& child : this->getChildren()) {
            for (auto v : child->getMembers()) {
                visited[v] = 1;
            }
        }
        bool contains_all_children = true;
        for (auto v : this->getMembers()) {
            if (!visited[v]) contains_all_children &= false;
        }
        return contains_all_children;
    }
};


class ITree {
public:
    virtual ~ITree() = default;

    // Root node of the decomposition
    virtual ITreeNode* getRoot() const = 0;
};


static bool same_members(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size()) return false;
    std::vector<int> aa = a, bb = b;
    std::sort(aa.begin(), aa.end());
    std::sort(bb.begin(), bb.end());
    return aa == bb;
}
#endif //OBLIVIOUSROUTING_RAECKE_TREE_H