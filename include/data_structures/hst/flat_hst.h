//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_FLAT_HST_H
#define OBLIVIOUSROUTING_FLAT_HST_H

#include <algorithm>
#include <vector>

struct HSTNodeFlat {
    int center;
    int members_begin, members_end;  // slice into FlatHST::members[]
    int children_begin, children_end; // slice into FlatHST::children_idx[]
};

struct FlatHST {
    std::vector<HSTNodeFlat> nodes;        // index 0 = root
    std::vector<int>         children_idx; // flat list of child node indices
    std::vector<int>         members;      // flat list of all member vertices

    // zero-copy child iteration: returns pointer pair into children_idx
    struct ChildRange { const int* b; const int* e;
        const int* begin() const { return b; }
        const int* end()   const { return e; }
    };
    ChildRange children(int i) const;

    // zero-copy member iteration
    struct MemberRange {
        const int* b; const int* e;
        const int* begin() const { return b; }
        const int* end()   const { return e; }
        int size()  const { return static_cast<int>(e - b); }
        bool empty() const { return b == e; }
        int operator[](int k) const { return b[k]; }
    };
    MemberRange memberRange(int i) const;
    int root()     const;
    int numNodes() const;
    bool isLeaf(int i) const;

    void print();
};


struct BuildNode {
    int center = -1;
    int parent = -1;           // index into build_nodes_, -1 = root
    std::vector<int> children;
    std::vector<int> members;
};

/** HSTBuilder: builds a mutable tree level-by-level from tree_oracle,
* then packs it into a flat FlatHST in one pass.
*/
class HSTBuilder {
public:

    explicit HSTBuilder(int n) {
        build_nodes_.reserve(4 * n);
        leaf_ids_.resize(n);
        for (int v = 0; v < n; ++v) {
            int idx = (int)build_nodes_.size();
            build_nodes_.push_back({v, -1, {}, {v}});
            leaf_ids_[v] = idx;
        }
    }

    int leafOf(int v) const;
    int size()        const;

    int addNode(int center);

    BuildNode& node(int i);
    const BuildNode& node(int i) const;

    /**
     * attach child to parent, propagate members upward
     */
    void attach(int parent_idx, int child_idx);

    void sortMembers(int i);

    /**
     * cleanUp + pack into flat FlatHST
     */
    FlatHST finalise(int root_idx);

    /**
     * if children and parent are the same, compress them into a single node
     */
    void compressTree();

private:
    std::vector<BuildNode> build_nodes_;
    std::vector<int>       leaf_ids_;

    void cleanUp(int i);
    static bool sameMembers(const std::vector<int>& a, const std::vector<int>& b);
};

int calculateFlatHSTHeight(const FlatHST& hst);

#endif //OBLIVIOUSROUTING_FLAT_HST_H