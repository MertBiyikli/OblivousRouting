//
// Created by Mert Biyikli on 17.09.25.
//

#ifndef OBLIVIOUSROUTING_MST_H
#define OBLIVIOUSROUTING_MST_H

#include <vector>
#include <tuple>
#include <random>
#include <algorithm>
#include <queue>
#include <unordered_map>
#include <limits>
#include "../datastructures/GraphCSR.h"
#include "hst.h"



/* Union find data structure to be used in the Kruskal's algorithm */
// ---------- DSU ----------
struct DSU {
    std::vector<int> p, r;
    explicit DSU(int n): p(n), r(n,0) {
        for (int i=0; i<n; ++i) p[i] = i;
    }
    int find(int x) { return p[x]==x ? x : p[x]=find(p[x]); }
    bool unite(int a, int b) {
        a = find(a); b = find(b);
        if (a==b) return false;
        if (r[a]<r[b]) std::swap(a,b);
        p[b]=a;
        if (r[a]==r[b]) ++r[a];
        return true;
    }
};

// ---------- MST tree struct ----------
struct MSTTree {
    int root;                                  // root vertex
    std::vector<int> parent;                   // parent[v] = parent of v (root points to itself)
    std::vector<std::vector<int>> children;    // adjacency list of MST as rooted tree
};

// ---------- Random MST oracle ----------
class RandomMST {
public:
    RandomMST(IGraph& g) : graph(g) {
        setGraph(g);
    }
    IGraph& graph;


    int n;
    std::vector<std::pair<int,int>> edges;
    std::vector<std::vector<int>> adj; // adjacency list of current MST
    std::vector<std::tuple<double,int,int>> keyed;

    std::vector<std::vector<double>> weights;
    void setGraph(const IGraph& g) {
        // clear all first
        edges.clear();
        keyed.clear();
        weights.clear();

        n = g.getNumNodes();

        for (int u = 0; u < n; ++u) {
            for (int v : g.neighbors(u)) {
                edges.emplace_back(u, v);
                keyed.emplace_back(g.getEdgeDistance(u, v), u, v);
            }
        }
    }


    void updateEdgeDistances(const std::vector<double>& distances) {
        keyed.clear();
        for (const auto& [u,v] : edges) {
            int e = graph.getEdgeId(u,v);
            double w = distances[e];
            keyed.emplace_back(w, u, v);
            this->graph.updateEdgeDistance(e, w);
        }
    }

    // Build a random MST edge set using Kruskal with random priorities
    std::vector<std::pair<int,int>> computeMST() {
        if (edges.empty() || keyed.empty()) return {};

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



/*
    // Turn MST edges into a rooted MSTTree (root at 0 by default)
    MSTTree buildMSTTree(const std::vector<std::pair<int,int>>& mst_edges,
                              int root=0) {
        adj = std::vector<std::vector<int>>(n);
        for (auto [u,v] : mst_edges) {
            adj[u].push_back(v);
            adj[v].push_back(u);
        }

        std::vector<int> parent(n, -1);
        std::vector<std::vector<int>> children(n);

        std::queue<int> q;
        q.push(root);
        parent[root] = root;

        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int v : adj[u]) {
                if (parent[v] != -1) continue;
                parent[v] = u;
                children[u].push_back(v);
                q.push(v);
            }
        }
        MSTTree t = {root, parent, children};
        return t;
    }
*/
    std::shared_ptr<HSTNode> buildHST(const std::vector<std::pair<int,int>>& mst_edges,
                              int root=0) {
        //auto mst_tree = buildMSTTree(mst_edges, root);
        std::vector<std::shared_ptr<HSTNode>> cluster(n);
        for (int v = 0; v < n; ++v) {
            cluster[v] = std::make_shared<HSTNode>(v);
            cluster[v]->id = v;
            cluster[v]->members = {v};
            cluster[v]->center = v;
        }

        // build Raecke Tree given the original MST tree
        DSU dsu(n);

        for (auto& [w, u, v] : keyed) {
            // if ( u >= v) continue; // process each edge once
            int pu = dsu.find(u);
            int pv = dsu.find(v);
            if (pu == pv) continue;

            auto parent = std::make_shared<HSTNode>(std::numeric_limits<int>::max());
            parent->children = {cluster[pu], cluster[pv]};
            cluster[pu]->parent = parent;
            cluster[pv]->parent = parent;

            // merge members
            parent->members = cluster[pu]->members;
            parent->members.insert(
                parent->members.end(),
                cluster[pv]->members.begin(),
                cluster[pv]->members.end()
            );

            dsu.unite(pu, pv);
            int new_rep = dsu.find(pu);
            parent->center = cluster[pu]->center; // arbitrary
            cluster[new_rep] = parent;
        }
        int rep = dsu.find(root);
        assert(cluster[dsu.find(root)] != nullptr);
        assert(cluster[dsu.find(root)]->children.size() > 0);
        return cluster[rep];
    }

};

#endif // OBLIVIOUSROUTING_MST_H
