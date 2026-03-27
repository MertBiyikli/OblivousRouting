#include "../../../../../../include/algorithms/mwu/oracle/tree/mst/mst_algo.h"
#include "../../../../../../include/data_structures/union_find/union_find.h"



MST::MST(IGraph& g) {
    keyed.clear();

    n = g.getNumNodes();

    for (int u = 0; u < n; ++u) {
        for (int v : g.neighbors(u)) {
            keyed.emplace_back(g.getEdgeDistance(u, v), u, v);
        }
    }
}

// Build a random MST edge set using Kruskal with random priorities
std::vector<std::pair<int,int>> MST::computeMST() {
    if ( keyed.empty()) return {};

    // apply kruskals algorithm
    std::sort(keyed.begin(), keyed.end(),
              [](auto& a, auto& b){ return std::get<0>(a) < std::get<0>(b); });


    DSU dsu(n);
    std::vector<std::pair<int,int>> mst;
    mst.reserve(n-1);

    for (auto& [key,u,v] : keyed) {
        if (dsu.unite(u,v)) {
            mst.emplace_back(u,v);
            if ((int)mst.size()+1 == n) break;
        }
    }
    return mst;
}


