//
// Created by Mert Biyikli on 24.03.26.
//

#ifndef OBLIVIOUSROUTING_UNION_FIND_H
#define OBLIVIOUSROUTING_UNION_FIND_H

#include <vector>
/* Union find data structure to be used in the Kruskal's algorithm */
class DSU {
public:
    std::vector<int> p, r;
    explicit DSU(int n): p(n), r(n,0) {
        for (int i=0; i<n; ++i) p[i] = i;
    }
    int find(int x);
    bool unite(int a, int b);
};

#endif //OBLIVIOUSROUTING_UNION_FIND_H