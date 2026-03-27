//
// Created by Mert Biyikli on 24.03.26.
//

#include "../../include/data_structures/union_find/union_find.h"
#include <algorithm>

int DSU::find(int x) {
    return p[x]==x ? x : p[x]=find(p[x]);
}

bool DSU::unite(int a, int b) {
    int a_rep = find(a);
    int b_rep = find(b);
    if (a_rep==b_rep) return false;
    if (r[a_rep]<r[b_rep]) {
        std::swap(a_rep,b_rep);
    }
    p[b_rep]=a_rep;
    if (r[a_rep]==r[b_rep]) {
        ++r[a_rep];
    }
    return true;
}