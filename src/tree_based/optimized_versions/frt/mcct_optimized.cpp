//
// Created by Mert Biyikli on 06.11.25.
//

#include "mcct_optimized.h"
#include "../../../utils/hash.h"

MCCT_optimized::MCCT_optimized() = default;

void MCCT_optimized::setGraph(const Graph_csr& g) {
    if (g.getNumNodes() <= 0 || g.getNumEdges() <= 0)
        throw std::runtime_error("MCCT_optimized::setGraph: empty graph");
    G = g;
}

void MCCT_optimized::setSeed(uint64_t seed) { rng.seed(seed); }
void MCCT_optimized::clearDemands() { demands.clear(); }


/*
double MCCT_optimized::estimateDiameter() const {
    // Lightweight 2-source Dijkstra heuristic on weighted graph.
    auto oneShot = [&](int src) {
        const int n = G.getNumNodes();
        std::vector<double> dist(n, std::numeric_limits<double>::infinity());
        using P = std::pair<double,int>;
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        dist[src]=0; pq.emplace(0,src);
        while(!pq.empty()){
            auto [d,u]=pq.top(); pq.pop();
            if (d!=dist[u]) continue;
            for (int eid =G.head[u]; eid < G.head[u+1]; ++eid) {
                int v = G.edges[eid].second;
                double w = G.getEdgeDistance(eid);
                if (dist[v] > d + w) { dist[v]=d+w; pq.emplace(dist[v], v); }
            }
        }
        int arg=-1; double far=-1;
        for (int i=0;i<n;++i) if (dist[i]<1e308 && dist[i]>far){far=dist[i]; arg=i;}
        return std::pair<int,double>(arg, far);
    };

    int start = 0;
    auto [mid, _] = oneShot(start);
    auto [end, diam] = oneShot(mid);
    if (!std::isfinite(diam) || diam<=0) diam = 1.0; // fallback
    return diam;
}
*/

/*
std::vector<int> MCCT_optimized::randomPermutation(int n) {
    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 0);
    std::shuffle(p.begin(), p.end(), rng);
    return p;
}

double MCCT_optimized::randomAlpha() {
    std::uniform_real_distribution<double> U(1.0, 2.0);
    return U(rng);
}
*/
void MCCT_optimized::truncatedDijkstra(int src, double R,
                                       std::vector<double>& dist,
                                       std::vector<int>& owner,
                                       int ownerId) const {
    using P = std::pair<double,int>;
    std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
    dist[src] = 0.0;
    pq.emplace(0.0, src);

    while(!pq.empty()){
        auto [d,u] = pq.top(); pq.pop();
        if (d!=dist[u]) continue;
        if (d > R) break;                 // beyond ball
        if (owner[u] != -1) continue;     // already claimed
        owner[u] = ownerId;

        for (int eid =G.head[u]; eid < G.head[u+1]; ++eid) {
            int v = G.edges[eid].second;
            double w = G.getEdgeDistance(eid);
            double nd = d + w;
            if (nd < dist[v] && nd <= R) {
                dist[v] = nd;
                pq.emplace(nd, v);
            }
        }
    }
}

int MCCT_optimized::makeLevel(double radius,
                              const std::vector<int>& perm,
                              std::vector<int>& owner,
                              std::vector<std::vector<int>>& clusters) const {
    const int n = G.getNumNodes();
    owner.assign(n, -1);
    int cid = 0;

    for (int pivot : perm) {
        if (owner[pivot] != -1) continue;

        // local distances (truncated)
        std::vector<double> dist(n, std::numeric_limits<double>::infinity());
        truncatedDijkstra(pivot, radius, dist, owner, cid);
        ++cid;
    }

    clusters.assign(cid, {});
    for (int v=0; v<n; ++v) {
        int c = owner[v];
        if (c>=0) clusters[c].push_back(v);
    }
    return cid;
}

void MCCT_optimized::ensureDemands() {
    if (!demands.empty()) return;
    int count = 0;
    for (int edge_id = 0; edge_id < G.getNumEdges(); ++edge_id) {
        const auto& [s, t] = G.edges[edge_id];
        demands.push_back({s, t, G.getEdgeCapacity(edge_id)});
        if (max_pairs > 0 && ++count >= max_pairs) return;
    }
}

void MCCT_optimized::buildTreeFromLevels(
    const std::vector<std::vector<std::vector<int>>>& levels,
    FRT_Tree& T) const
{
    // Assumes levels[ℓ] is vector of clusters (each cluster is vector<int> of vertices),
    // ℓ = L..0 (top..bottom). Create a tree node per cluster and connect to its parent.

    const int L = (int)levels.size();
    if (L==0) return;

    // Map each vertex to its cluster id at each level to find parents quickly.
    const int n = G.getNumNodes();
    std::vector<std::vector<int>> who(n, std::vector<int>(L, -1)); // who[v][ℓ] = cluster id at level ℓ

    for (int ell = 0; ell < L; ++ell) {
        for (int cid = 0; cid < (int)levels[ell].size(); ++cid) {
            for (int v : levels[ell][cid]) who[v][ell] = cid;
        }
    }

    // Create nodes per cluster per level
    std::vector<std::vector<int>> nodeId(L);
    for (int ell=0; ell<L; ++ell) {
        nodeId[ell].resize(levels[ell].size(), -1);
        for (int cid=0; cid<(int)levels[ell].size(); ++cid) {
            nodeId[ell][cid] = T.addClusterNode(levels[ell][cid], ell); // You have this in your FRT_Tree
        }
    }

    // Connect child level to parent level
    for (int ell=0; ell<L-1; ++ell) {
        for (int cid=0; cid<(int)levels[ell].size(); ++cid) {
            // Pick any vertex in this cluster to find parent cluster at level ell+1
            int v = levels[ell][cid].empty() ? -1 : levels[ell][cid][0];
            if (v < 0) continue;
            int parent_cid = who[v][ell+1];
            if (parent_cid >= 0) {
                T.addTreeEdge(nodeId[ell][cid], nodeId[ell+1][parent_cid]);
            }
        }
    }

    // Optionally, you can set edge lengths per level if your FRT_Tree supports scale
    // e.g., T.setLevelLength(ell, 2^ell) – depends on your Tree API.
}

FRT_Tree MCCT_optimized::build(bool debug) {
    if (G.getNumNodes() <= 0) throw std::runtime_error("build(): graph not set");
    ensureDemands();

    const int n = G.getNumNodes();
    double diam = estimateDiameter();
    // L such that 2^L ≳ diam
    int L = (int)std::ceil(std::log2(std::max(1.0, diam))) + 1;

    auto perm  = randomPermutation(n);
    double a   = randomAlpha();           // α in [1,2)

    // Build levels from top (L-1) down to 0
    std::vector<std::vector<std::vector<int>>> levels; // levels[ell] -> clusters, ell=0..L-1 (bottom..top)
    levels.resize(L);

    // We’ll fill top -> down, and later reverse to bottom..top for tree builder
    std::vector<int> owner;
    for (int ell = L-1; ell >= 0; --ell) {
        double radius = a * std::ldexp(1.0, ell);   // α·2^ell
        std::vector<std::vector<int>> clusters;
        makeLevel(radius, perm, owner, clusters);
        levels[ell] = std::move(clusters);
        if (debug) {
            std::cout << "Level " << ell << " radius=" << radius
                      << " clusters=" << levels[ell].size() << "\n";
        }
    }

    FRT_Tree T;
    buildTreeFromLevels(levels, T);
    return T;
}