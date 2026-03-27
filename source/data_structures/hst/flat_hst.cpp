//
// Created by Mert Biyikli on 20.03.26.
//

#include "../../include/data_structures/hst/flat_hst.h"
#include <unordered_map>

FlatHST::ChildRange FlatHST::children(int i) const {
    const auto& n = nodes[i];
    if (n.children_begin < 0) return {nullptr, nullptr};
    return { children_idx.data() + n.children_begin,
             children_idx.data() + n.children_end };
}

FlatHST::MemberRange FlatHST::memberRange(int i) const {
    const auto& n = nodes[i];
    return { members.data() + n.members_begin,
             members.data() + n.members_end };
}

int FlatHST::root() const {
    return 0;
}

bool FlatHST::isLeaf(int i) const {
    return nodes[i].children_begin == nodes[i].children_end;
}

int FlatHST::numNodes() const {
    return static_cast<int>(nodes.size());
}

void FlatHST::print() {
    for (int i = 0; i < numNodes(); ++i) {
        const auto& n = nodes[i];
        printf("Node %d: center=%d, members=[", i, n.center);
        for (int j = n.members_begin; j < n.members_end; ++j) {
            printf("%d ", members[j]);
        }
        printf("], children=[");
        for (int j = n.children_begin; j < n.children_end; ++j) {
            printf("%d ", children_idx[j]);
        }
        printf("]\n");
    }
}


int HSTBuilder::leafOf(int v) const {
    return leaf_ids_[v];
}

int HSTBuilder::size() const {
    return static_cast<int>(build_nodes_.size());
}

int HSTBuilder::addNode(int center) {
    int idx = (int)build_nodes_.size();
    build_nodes_.push_back({center, -1, {}, {}});
    return idx;
}

BuildNode& HSTBuilder::node(int i) {
    return build_nodes_[i];
}

const BuildNode& HSTBuilder::node(int i) const {
    return build_nodes_[i];
}


void HSTBuilder::attach(int parent_idx, int child_idx) {
    BuildNode& parent = build_nodes_[parent_idx];
    BuildNode& child  = build_nodes_[child_idx];
    parent.children.push_back(child_idx);
    child.parent = parent_idx;
    // propagate members upward
    parent.members.insert(parent.members.end(), child.members.begin(), child.members.end());
}


void HSTBuilder::sortMembers(int i) {
    auto& m = build_nodes_[i].members;
    std::sort(m.begin(), m.end());
    m.erase(std::unique(m.begin(), m.end()), m.end());
}

void HSTBuilder::compressTree() {
    for (int i = 0; i < size(); ++i) {
        BuildNode& nd = build_nodes_[i];
        if (nd.members.size() <= 1) continue;
        for (int c : nd.children) {
            BuildNode& cn = build_nodes_[c];
            if (sameMembers(nd.members, cn.members)) {
                // absorb child: adopt grandchildren
                for (int gc : cn.children) {
                    build_nodes_[gc].parent = i;
                    nd.children.push_back(gc);
                }
                cn.children.clear();
            }
        }
    }
    // remove absorbed children
    for (int i = 0; i < size(); ++i) {
        BuildNode& nd = build_nodes_[i];
        std::vector<int> new_children;
        new_children.reserve(nd.children.size());
        for (int c : nd.children) {
            if (!build_nodes_[c].children.empty()) new_children.push_back(c);
        }
        nd.children.swap(new_children);
    }
}

bool HSTBuilder::sameMembers(const std::vector<int> &a, const std::vector<int> &b) {
    if (a.size() != b.size()) return false;
    std::vector<int> aa = a, bb = b;
    std::sort(aa.begin(), aa.end());
    std::sort(bb.begin(), bb.end());
    return aa == bb;
}


void HSTBuilder::cleanUp(int i) {
    BuildNode& nd = build_nodes_[i];
    if (nd.members.empty()) return;
    if (nd.members.size() == 1) { nd.center = nd.members[0]; nd.children.clear(); return; }

    std::vector<int> kids = nd.children;
    for (int c : kids) cleanUp(c);

    std::vector<int> new_children;
    new_children.reserve(nd.children.size());
    for (int c : nd.children) {
        BuildNode& cn = build_nodes_[c];
        if (sameMembers(nd.members, cn.members)) {
            for (int gc : cn.children) { build_nodes_[gc].parent = i; new_children.push_back(gc); }
            cn.children.clear();
        } else {
            new_children.push_back(c);
        }
    }
    nd.children.swap(new_children);
}

// Inline implementation of finalise
FlatHST HSTBuilder::finalise(int root_idx) {
    cleanUp(root_idx);

    // BFS to assign new contiguous indices (root = 0)
    std::vector<int> old_to_new(build_nodes_.size(), -1);
    std::vector<int> bfs;
    bfs.reserve(build_nodes_.size());
    bfs.push_back(root_idx);
    old_to_new[root_idx] = 0;
    for (int qi = 0; qi < (int)bfs.size(); ++qi) {
        for (int c : build_nodes_[bfs[qi]].children) {
            old_to_new[c] = (int)bfs.size();
            bfs.push_back(c);
        }
    }

    FlatHST hst;
    const int N = (int)bfs.size();
    hst.nodes.resize(N);

    // pack members
    for (int ni = 0; ni < N; ++ni) {
        const BuildNode& bn = build_nodes_[bfs[ni]];
        auto& hn = hst.nodes[ni];
        hn.center       = bn.center;
        hn.members_begin = (int)hst.members.size();
        hst.members.insert(hst.members.end(), bn.members.begin(), bn.members.end());
        hn.members_end   = (int)hst.members.size();
    }

    // pack children
    for (int ni = 0; ni < N; ++ni) {
        const BuildNode& bn = build_nodes_[bfs[ni]];
        auto& hn = hst.nodes[ni];
        if (bn.children.empty()) {
            hn.children_begin = hn.children_end = 0;
        } else {
            hn.children_begin = (int)hst.children_idx.size();
            for (int c : bn.children) hst.children_idx.push_back(old_to_new[c]);
            hn.children_end = (int)hst.children_idx.size();
        }
    }
    return hst;
}


// Helper function to calculate height of flat HST
int calculateFlatHSTHeight(const FlatHST& hst) {
    if (hst.nodes.empty()) return 0;

    // Start from the root and traverse the tree structure
    int maxHeight = 0;
    std::vector<int> stack;
    stack.push_back(hst.root());

    std::unordered_map<int, int> heights;

    // Post-order traversal to calculate heights
    std::vector<bool> visited(hst.nodes.size(), false);
    std::vector<int> path;

    while (!stack.empty()) {
        int idx = stack.back();

        if (visited[idx]) {
            path.pop_back();
            stack.pop_back();

            int h = 1;
            auto range = hst.children(idx);
            for (auto it = range.begin(); it != range.end(); ++it) {
                int child_idx = *it;
                h = std::max(h, 1 + heights[child_idx]);
            }
            heights[idx] = h;
            maxHeight = std::max(maxHeight, h);
        } else {
            visited[idx] = true;
            path.push_back(idx);
            auto range = hst.children(idx);
            for (auto it = range.begin(); it != range.end(); ++it) {
                int child_idx = *it;
                if (!visited[child_idx]) {
                    stack.push_back(child_idx);
                }
            }
        }
    }

    return maxHeight;
}
