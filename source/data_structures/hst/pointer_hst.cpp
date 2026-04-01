//
// Created by Mert Biyikli on 20.03.26.
//

#include "../../include/data_structures/hst/pointer_hst.h"
#include <algorithm>

bool same_members(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size()) return false;
    std::vector<int> aa = a, bb = b;
    std::sort(aa.begin(), aa.end());
    std::sort(bb.begin(), bb.end());
    return aa == bb;
}

// Helper function to calculate height of pointer HST
int calculatePointerHSTHeight(const std::shared_ptr<HSTNode>& node) {
    if (!node) return 0;
    if (node->children.empty()) return 1;
    int maxChildHeight = 0;
    for (const auto& child : node->children) {
        maxChildHeight = std::max(maxChildHeight, calculatePointerHSTHeight(child));
    }
    return 1 + maxChildHeight;
}