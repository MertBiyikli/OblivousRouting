//
// Created by Mert Biyikli on 07.03.26.
//

#ifndef OBLIVIOUSROUTING_FLAT_HST_H
#define OBLIVIOUSROUTING_FLAT_HST_H

#include <vector>

struct HSTNodeFlat {
    int center;
    int members_begin, members_end;  // slice into HST::members[]
    int children_begin, children_end; // slice into HST::children_idx[]
};

struct HST {
    std::vector<HSTNodeFlat> nodes;        // index 0 = root
    std::vector<int>         children_idx; // flat list of child node indices
    std::vector<int>         members;      // flat list of all member vertices

    // zero-copy child iteration: returns pointer pair into children_idx
    struct ChildRange { const int* b; const int* e;
        const int* begin() const { return b; }
        const int* end()   const { return e; }
    };
    ChildRange children(int i) const {
        const auto& n = nodes[i];
        if (n.children_begin < 0) return {nullptr, nullptr};
        return { children_idx.data() + n.children_begin,
                 children_idx.data() + n.children_end };
    }

    // zero-copy member iteration
    struct MemberRange { const int* b; const int* e;
        const int* begin() const { return b; }
        const int* end()   const { return e; }
        int size()  const { return static_cast<int>(e - b); }
        bool empty() const { return b == e; }
        int operator[](int k) const { return b[k]; }
    };
    MemberRange memberRange(int i) const {
        const auto& n = nodes[i];
        return { members.data() + n.members_begin,
                 members.data() + n.members_end };
    }

    int root()     const { return 0; }
    int numNodes() const { return static_cast<int>(nodes.size()); }
    bool isLeaf(int i) const {
        return nodes[i].children_begin == nodes[i].children_end;
    }
};

// HSTBuilder: builds a mutable tree level-by-level from tree_oracle,
// then packs it into a flat HST in one pass.
class HSTBuilder {
public:
    struct BuildNode {
        int center = -1;
        int parent = -1;           // index into build_nodes_, -1 = root
        std::vector<int> children;
        std::vector<int> members;
    };

    explicit HSTBuilder(int n) {
        build_nodes_.reserve(4 * n);
        leaf_ids_.resize(n);
        for (int v = 0; v < n; ++v) {
            int idx = (int)build_nodes_.size();
            build_nodes_.push_back({v, -1, {}, {v}});
            leaf_ids_[v] = idx;
        }
    }

    int leafOf(int v) const { return leaf_ids_[v]; }
    int size()        const { return (int)build_nodes_.size(); }

    int addNode(int center) {
        int idx = (int)build_nodes_.size();
        build_nodes_.push_back({center, -1, {}, {}});
        return idx;
    }

    BuildNode& node(int i)       { return build_nodes_[i]; }
    const BuildNode& node(int i) const { return build_nodes_[i]; }

    // attach child to parent, propagate members upward
    void attach(int parent_idx, int child_idx) {
        build_nodes_[child_idx].parent = parent_idx;
        build_nodes_[parent_idx].children.push_back(child_idx);
        auto& pm = build_nodes_[parent_idx].members;
        auto& cm = build_nodes_[child_idx].members;
        pm.insert(pm.end(), cm.begin(), cm.end());
    }

    void sortMembers(int i) {
        auto& m = build_nodes_[i].members;
        std::sort(m.begin(), m.end());
        m.erase(std::unique(m.begin(), m.end()), m.end());
    }

    // cleanUp + pack into flat HST
    HST finalise(int root_idx);

private:
    std::vector<BuildNode> build_nodes_;
    std::vector<int>       leaf_ids_;

    void cleanUp(int i);
    static bool sameMembers(const std::vector<int>& a, const std::vector<int>& b) {
        return a.size() == b.size() && a == b; // both pre-sorted
    }
};

// Implementation of finalise() and cleanUp() inline:

inline void HSTBuilder::cleanUp(int i) {
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

inline HST HSTBuilder::finalise(int root_idx) {
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

    HST hst;
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


#endif //OBLIVIOUSROUTING_FLAT_HST_H