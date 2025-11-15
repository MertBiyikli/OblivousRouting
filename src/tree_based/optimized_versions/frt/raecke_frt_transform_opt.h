//
// Created by Mert Biyikli on 05.11.25.
//

#ifndef OBLIVIOUSROUTING_RAECKE_FRT_TRANSFORM_OPT_H
#define OBLIVIOUSROUTING_RAECKE_FRT_TRANSFORM_OPT_H


template<typename Tree>
class RaeckeFRTTransformOptimized{
public:

    void addTree(
            Tree& t,
            double           lambda,
            Graph_csr&           graph)
    {
        // For now, just a placeholder
        // Actual implementation would go here

    }

    std::unordered_map<
    std::pair<int,int>,
    std::unordered_map<std::pair<int,int>, double>
    > getRoutingTable() {
        // Placeholder
        return {};
    }
};

struct FlowStore {
    int m = 0;
    std::vector<uint64_t> st_idx;
    std::vector<double> val;
    std::vector<uint32_t> edge_ptr;
    std::vector<std::vector<std::pair<uint64_t, double>>> buckets;

    static inline uint64_t pack(uint32_t s, uint32_t t) {
        return (uint64_t(s) << 32) | uint64_t(t);
    }

    static inline std::pair<uint32_t, uint32_t> unpack(uint64_t id) {
        return {uint32_t(id >> 32), uint32_t(id & 0xffffffff)};
    }

    void begin(int M) {
        m = M;
        buckets.assign(m, {});
    }

    inline void add(int e, uint32_t s, uint32_t t, double x) {
        buckets[e].emplace_back(pack(s, t), x);
    }

    void finalize() {
        edge_ptr.assign(m + 1, 0);
        size_t total = 0;
        for (int e = 0; e < m; ++e) {
            auto &B = buckets[e];
            if (B.empty()) continue;
            std::sort(B.begin(), B.end(), [](auto &a, auto &b) { return a.first < b.first; });

            size_t w = 0;
            for (size_t i = 0; i < B.size();) {
                uint64_t k = B[i].first;
                double sumv = 0.0;
                size_t j = i;
                while (j < B.size() && B[j].first == k) { sumv += B[j].second; ++j; }
                B[w++] = {k, sumv};
                i = j;
            }
            B.resize(w);
            total += w;
            edge_ptr[e + 1] = total;
        }
        st_idx.resize(total);
        val.resize(total);

        size_t off = 0;
        for (int e = 0; e < m; ++e) {
            for (auto &p : buckets[e]) {
                st_idx[off] = p.first;
                val[off] = p.second;
                ++off;
            }
        }
        buckets.clear();
    }

    struct EdgeFlowRange {
        const FlowStore *fs;
        int e;
        struct It {
            const FlowStore *fs;
            size_t i;
            bool operator!=(const It &o) const { return i != o.i; }
            void operator++() { ++i; }
            auto operator*() const {
                return std::pair<uint64_t, double>(fs->st_idx[i], fs->val[i]);
            }
        };
        It begin() const { return {fs, fs->edge_ptr[e]}; }
        It end() const { return {fs, fs->edge_ptr[e + 1]}; }
    };

    EdgeFlowRange flows_on_edge(int e) const { return {this, e}; }
};

struct TreeAgg {
    int n = 0;
    std::vector<int> parent, depth, heavy, head, pos;
    std::vector<int> size;
    int cur = 0;

    struct BIT {
        int n;
        std::vector<double> bit;
        BIT(int n = 0) : n(n), bit(n + 2, 0.0) {}
        void add(int i, double v) { for (++i; i <= n; i += i & -i) bit[i] += v; }
        double sum(int i) const { double s = 0; for (++i; i > 0; i -= i & -i) s += bit[i]; return s; }
        void range_add(int l, int r, double v) { add(l, v); add(r + 1, -v); }
        double point(int i) const { return sum(i); }
    } bit;

    void build(const std::vector<std::vector<int>> &adj, int root = 0) {
        n = (int)adj.size();
        parent.assign(n, -1);
        depth.assign(n, 0);
        size.assign(n, 0);
        heavy.assign(n, -1);
        head.assign(n, 0);
        pos.assign(n, 0);

        std::function<int(int)> dfs = [&](int u) {
            int sz = 1, max_sub = 0;
            for (int v : adj[u]) if (v != parent[u]) {
                parent[v] = u;
                depth[v] = depth[u] + 1;
                int sub = dfs(v);
                if (sub > max_sub) { max_sub = sub; heavy[u] = v; }
                sz += sub;
            }
            return size[u] = sz;
        };
        parent[root] = -1;
        dfs(root);

        cur = 0;
        std::function<void(int,int)> decomp = [&](int u, int h) {
            head[u] = h;
            pos[u] = cur++;
            if (heavy[u] != -1) decomp(heavy[u], h);
            for (int v : adj[u]) if (v != parent[u] && v != heavy[u])
                decomp(v, v);
        };
        decomp(root, root);
        bit = BIT(n);
    }

    void add_path(int u, int v, double a) {
        while (head[u] != head[v]) {
            if (depth[head[u]] < depth[head[v]]) std::swap(u, v);
            int h = head[u];
            bit.range_add(pos[h], pos[u], a);
            u = parent[h];
        }
        if (depth[u] > depth[v]) std::swap(u, v);
        bit.range_add(pos[u] + 1, pos[v], a); // exclude LCA node itself
    }

    double edge_flow(int x) const { return bit.point(pos[x]); }
};

struct SchemeApplier {
    struct WeightedTree {
        double lambda;
        TreeAgg tree;
        std::vector<std::vector<int>> adj;
    };
    std::vector<WeightedTree> trees;

    void addTree(double lambda, const std::vector<std::vector<int>>& adj, int root = 0) {
        trees.push_back({});
        trees.back().lambda = lambda;
        trees.back().adj = adj;
        trees.back().tree.build(adj, root);
    }

    // Apply all demands to all trees, return total edge load per original node index
    std::vector<double> apply(const std::vector<std::tuple<int,int,double>>& demands) {
        if (trees.empty()) return {};
        int n = trees[0].tree.n;
        std::vector<double> agg(n, 0.0);

        for (auto& wT : trees) {
            auto& T = wT.tree;
            for (auto [s, t, d] : demands)
                T.add_path(s, t, d);
            for (int v = 1; v < n; ++v) {
                double f = T.edge_flow(v);
                agg[v] += f * wT.lambda;
            }
        }
        return agg;
    }
};

#endif //OBLIVIOUSROUTING_RAECKE_FRT_TRANSFORM_OPT_H