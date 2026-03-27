
#ifndef OBLIVIOUSROUTING_TREE_ORACLE_H
#define OBLIVIOUSROUTING_TREE_ORACLE_H

#include <list>

#include "hst.h"
#include "flat_hst.h"
#include "../datastructures/IGraph.h"
#include "utils/ultrametric_tree.h"
#include "utils/quotient_graph.h"
#include "mst.h"
#include <vector>

// ---------------------------------------------------------------------------
// TreeOracle<T>
//
// T = std::shared_ptr<HSTNode>  → pointer-based tree (computeHST path)
// T = FlatHST                       → flat array tree   (getFlatTree path)
//
// Subclasses (FRT, FastCKR, RandomMST) only override computeLevelPartition,
// which is independent of T, so they become TreeOracle<T> subclasses with
// no changes to their own code.
// ---------------------------------------------------------------------------
template<typename T>
class TreeOracle {
public:
    explicit TreeOracle(IGraph& graph) : graph(graph) {
        n = graph.getNumNodes();
        diameter = graph.GetDiameter();
        applyMendelScaling = false;
    }

    TreeOracle(IGraph& graph, bool activateMendelScaling) : graph(graph) {
        n = graph.getNumNodes();
        diameter = graph.GetDiameter();
        this->applyMendelScaling = activateMendelScaling;
        total_time_spent_on_mendel_scaling = 0.0;
    }
    virtual ~TreeOracle() = default;

    IGraph& graph;
    int n;
    double diameter;
    std::vector<double> scales;
    std::vector<int> perm;

    bool applyMendelScaling = false;
    MendelScaling::UltrametricTree ultrametric;
    MendelScaling::QuotientConstruction qc;
    double total_time_spent_on_mendel_scaling;

    // Single entry point — dispatches at compile time based on T.
    virtual T getTree(std::vector<double>& distances) {
        total_time_spent_on_mendel_scaling = 0;
        updateDistances(distances);
        if constexpr (std::is_same_v<T, std::shared_ptr<HSTNode>>) {
            return computeHST();
        } else {
            return computeFlatHST();
        }
    }

    // -----------------------------------------------------------------------
    // Pointer-based HST
    // -----------------------------------------------------------------------
    std::shared_ptr<HSTNode> computeHST() {
        preprocess();

        std::vector<std::shared_ptr<HSTNode>> prev_nodes;
        for (int v = 0; v < n; ++v) {
            prev_nodes.emplace_back(std::make_shared<HSTNode>(v));
            prev_nodes.back()->id = v;
            prev_nodes.back()->members = {v};
            prev_nodes.back()->center = v;
        }

        computeScales();

        std::shared_ptr<HSTNode> root = std::make_shared<HSTNode>(std::numeric_limits<int>::max());

        if (applyMendelScaling) {
            auto start = timeNow();
            qc.preprocessEdges(graph);
            total_time_spent_on_mendel_scaling += duration(timeNow() - start);
        }

        std::vector<std::shared_ptr<HSTNode>> current_tree_nodes = prev_nodes;

        for (const double Delta : scales) {
            std::vector<std::shared_ptr<HSTNode>> current_level;
            HSTLevel L;

            MendelScaling::QuotientLevel Q;
            if (applyMendelScaling) {
                auto start = timeNow();
                Q = qc.constructQuotientGraph(ultrametric, Delta, graph);
                if (Q.Gq->getNumNodes() <= 1) continue;
                total_time_spent_on_mendel_scaling += duration(timeNow() - start);
                computeQuotientLevelPartition(Q, L, Delta);
            } else {
                computeLevelPartition(graph, L, perm, Delta);
            }
            current_level = buildTreeLevel(prev_nodes, L);

            for (auto& node : current_level)
                current_tree_nodes.push_back(node);

            prev_nodes = current_level;
        }
        finishTree(root, current_tree_nodes);
        return root;
    }

    // -----------------------------------------------------------------------
    // Flat FlatHST
    // -----------------------------------------------------------------------
    FlatHST computeFlatHST() {
        preprocess();

        HSTBuilder builder(n);
        std::vector<int> prev_level(n);
        for (int v = 0; v < n; ++v)
            prev_level[v] = builder.leafOf(v);

        if (applyMendelScaling) {
            auto start = timeNow();
            qc.preprocessEdges(graph);
            total_time_spent_on_mendel_scaling += duration(timeNow() - start);
        }

        computeScales();
        for (const double Delta : scales) {
            HSTLevel L;
            if (applyMendelScaling) {
                MendelScaling::QuotientLevel Q;
                auto start = timeNow();
                Q = qc.constructQuotientGraph(ultrametric, Delta, graph);
                if (Q.Gq->getNumNodes() <= 1) continue;
                total_time_spent_on_mendel_scaling += duration(timeNow() - start);
                computeQuotientLevelPartition(Q, L, Delta);
            } else {
                computeLevelPartition(graph, L, perm, Delta);
            }
            prev_level = buildFlatTreeLevel(builder, prev_level, L);
        }

        int root_idx = builder.addNode(builder.node(prev_level[0]).center);
        std::vector<bool> attached(builder.size(), false);
        for (int v = 0; v < n; ++v) {
            int top = prev_level[v];
            if (attached[top]) continue;
            attached[top] = true;
            builder.attach(root_idx, top);
        }
        builder.sortMembers(root_idx);

        return builder.finalise(root_idx);
    }

    // -----------------------------------------------------------------------
    // The only pure virtual — subclasses only need to implement this.
    // -----------------------------------------------------------------------
    virtual void computeLevelPartition(IGraph& g, HSTLevel& level,
                                       const std::vector<int>& x_perm,
                                       double delta) = 0;

    // -----------------------------------------------------------------------
    // Shared helpers (unchanged from original)
    // -----------------------------------------------------------------------
    void preprocess() {
        n = graph.getNumNodes();
        if (n == 0) throw std::invalid_argument("The graph has no nodes.");
        diameter = graph.GetDiameter();

        perm.clear();
        for (int v = 0; v < n; ++v) perm.push_back(v);

        if (applyMendelScaling) {
            auto start = timeNow();
            RandomMST mst_oracle(graph);
            auto mst = mst_oracle.computeMST();
            std::vector<double> mst_weights;
            mst_weights.reserve(mst.size());
            for (auto [u, v] : mst)
                mst_weights.push_back(graph.getEdgeDistance(u, v));
            ultrametric.buildFromMST(graph.getNumNodes(), mst, mst_weights);
            total_time_spent_on_mendel_scaling += duration(timeNow() - start);
            assert(ultrametric.root != -1);
        }
    }

    static std::vector<int> buildRepresentatives(const std::vector<int>& sigma, int nQ) {
        const int nG = (int)sigma.size();
        std::vector<int> rep(nQ, -1);
        for (int v = 0; v < nG; ++v) {
            int q = sigma[v];
            if (q < 0 || q >= nQ)
                throw std::runtime_error("sigma[v] out of range at v=" + std::to_string(v));
            if (rep[q] == -1) rep[q] = v;
        }
        for (int q = 0; q < nQ; ++q)
            if (rep[q] == -1)
                throw std::runtime_error("Quotient node q=" + std::to_string(q) + " has no preimage");
        return rep;
    }

    void computeQuotientLevelPartition(MendelScaling::QuotientLevel& Q, HSTLevel& level, double delta) {
        HSTLevel qL;
        std::vector<int> x_perm;
        std::vector<char> used(Q.Gq->getNumNodes(), 0);

        for (int v : perm) {
            int vq = Q.sigma_compact_of_v[v];
            if (vq < 0 || vq >= Q.Gq->getNumNodes())
                throw std::runtime_error(
                    "computeQuotientLevelPartition: sigma_compact_of_v[v] out of range at v="
                    + std::to_string(v));
            if (!used[vq]) {
                x_perm.push_back(vq);
                used[vq] = 1;
            }
        }

        computeLevelPartition(*Q.Gq, qL, x_perm, delta / 2);

        const int nG = graph.getNumNodes();
        level.R = qL.R;
        level.centers.clear();
        level.owner.resize(nG);
        for (int v = 0; v < nG; ++v) level.owner[v] = v;

        // Deterministic quotient representatives induced by global permutation.
        std::vector<int> q_rep(Q.Gq->getNumNodes(), -1);
        for (int v : perm) {
            int q = Q.sigma_compact_of_v[v];
            if (q < 0 || q >= Q.Gq->getNumNodes())
                throw std::runtime_error(
                    "computeQuotientLevelPartition: sigma_compact_of_v[v] out of range at v="
                    + std::to_string(v));
            if (q_rep[q] == -1) q_rep[q] = v;
        }
        for (int q = 0; q < Q.Gq->getNumNodes(); ++q) {
            if (q_rep[q] == -1)
                throw std::runtime_error(
                    "computeQuotientLevelPartition: quotient node " + std::to_string(q) +
                    " has no representative in perm");
        }

        std::unordered_map<int, int> q_center_to_g_center;
        q_center_to_g_center.reserve(qL.centers.size());

        for (int q_center : qL.centers) {
            if (q_center < 0 || q_center >= Q.Gq->getNumNodes())
                throw std::runtime_error(
                    "computeQuotientLevelPartition: quotient center out of range: " +
                    std::to_string(q_center));

            int g_center = q_rep[q_center];
            level.centers.push_back(g_center);
            q_center_to_g_center[q_center] = g_center;
        }

        for (int vq = 0; vq < Q.Gq->getNumNodes(); ++vq) {
            int cluster_q = qL.owner[vq];
            if (cluster_q == -1)
                throw std::runtime_error(
                    "computeQuotientLevelPartition: node " + std::to_string(vq) + " has no cluster!");

            auto it = q_center_to_g_center.find(cluster_q);
            if (it == q_center_to_g_center.end())
                throw std::runtime_error(
                    "computeQuotientLevelPartition: cluster center " + std::to_string(cluster_q) +
                    " missing lifted representative");

            int cluster_g = it->second;
            for (int v : Q.members_of_q[vq]) level.owner[v] = cluster_g;
        }
    }


    void updateDistances(std::vector<double>& distances) {
        for (int e = 0; e < graph.getNumDirectedEdges(); ++e)
            graph.updateEdgeDistance(e, distances[e]);
    }

    void computeScales() {
        scales.clear();
        if (n == 0) throw std::invalid_argument("TreeOracle: The graph has no nodes.");
        if (applyMendelScaling) computeMendelScales();
        else                    computeNaiveScales();
    }

    std::vector<std::shared_ptr<HSTNode>> buildTreeLevel(
        std::vector<std::shared_ptr<HSTNode>>& prev_nodes,
        const HSTLevel& L) {
        std::vector<std::shared_ptr<HSTNode>> parents(L.centers.size());

        // center value -> index in parents
        std::unordered_map<int, int> center_to_idx;
        center_to_idx.reserve(L.centers.size());

        for (size_t c = 0; c < L.centers.size(); ++c) {
            parents[c] = std::make_shared<HSTNode>(-1);
            parents[c]->center = L.centers[c];
            parents[c]->parent.reset();
            parents[c]->children.clear();
            center_to_idx[L.centers[c]] = static_cast<int>(c);
        }

        for (const auto& child : prev_nodes) {
            if (!child) continue;

            int child_center = child->center;
            if (child_center < 0 || child_center >= static_cast<int>(L.owner.size())) {
                throw std::runtime_error(
                    "buildTreeLevel: child center out of range in L.owner: " +
                    std::to_string(child_center));
            }

            int cluster = L.owner[child_center];
            auto it = center_to_idx.find(cluster);
            if (it == center_to_idx.end()) {
                throw std::runtime_error(
                    "buildTreeLevel: cluster center not found in L.centers: " +
                    std::to_string(cluster));
            }

            auto& parent = parents[it->second];
            child->parent = parent;
            parent->children.push_back(child);
            parent->members.insert(parent->members.end(),
                                   child->members.begin(),
                                   child->members.end());
        }

        for (auto& p : parents) {
            auto& M = p->members;
            std::sort(M.begin(), M.end());
            M.erase(std::unique(M.begin(), M.end()), M.end());
        }

        std::vector<std::shared_ptr<HSTNode>> cleaned;
        cleaned.reserve(parents.size());
        for (auto& p : parents) {
            if (!p->members.empty()) cleaned.push_back(p);
        }

        return cleaned;
    }

    void computeMendelScales() {
        double min_distance = std::numeric_limits<double>::max();
        for (int e = 0; e < graph.getNumDirectedEdges(); e++) {
            double w = graph.getEdgeDistance(e);
            if (w < min_distance) min_distance = w;
        }
        if (diameter == 0)
            throw std::runtime_error("TreeOracle: Graph has zero diameter; cannot build tree.");
        double Delta = 1.0;
        while (Delta < diameter) Delta *= 8.0;
        for (; Delta >= min_distance; Delta /= 8.0) scales.push_back(Delta);
        std::reverse(scales.begin(), scales.end());
    }

    void computeNaiveScales() {
        if (diameter == 0)
            throw std::runtime_error("TreeOracle: Graph has zero diameter; cannot build tree.");
        int i = static_cast<int>(std::ceil(std::log2(diameter) / std::log2(2.0))) + 1;
        for (; i >= 0; --i) {
            double scale = (double)(1ull << i);
            if (scale <= diameter) scales.push_back(scale);
        }
        std::reverse(scales.begin(), scales.end());
    }

    void finishTree(std::shared_ptr<HSTNode>& root, std::vector<std::shared_ptr<HSTNode>>& prev_nodes) {
        if (prev_nodes.empty()) {
            root = std::make_shared<HSTNode>(std::numeric_limits<int>::max());
            root->members.resize(graph.getNumNodes());
            std::iota(root->members.begin(), root->members.end(), 0);
            root->center = (!root->members.empty() ? root->members[0] : -1);
        } else {
            root = std::make_shared<HSTNode>(std::numeric_limits<int>::max());
            root->members.resize(graph.getNumNodes());
            std::iota(root->members.begin(), root->members.end(), 0);
            root->center = prev_nodes.back()->center;

            std::set<int> members_set;
            for (int i = (int)prev_nodes.size() - 1; i >= 0; --i) {
                auto& node = prev_nodes[i];
                for (auto& v : node->members) members_set.insert(v);
                root->children.push_back(node);
                if (members_set.size() == (size_t)graph.getNumNodes()) break;
            }
        }
        cleanUpTree(root);
    }

    void cleanUpTree(std::shared_ptr<HSTNode>& node) {
        if (!node) return;
        if (node->members.empty()) {
            auto parent = node->parent.lock();
            if (parent) {
                auto& siblings = parent->children;
                siblings.erase(std::remove(siblings.begin(), siblings.end(), node), siblings.end());
            }
            return;
        }
        if (node->members.size() == 1) {
            node->center = node->members[0];
            node->children.clear();
        }
        for (auto& ch : node->children) cleanUpTree(ch);

        std::vector<std::shared_ptr<HSTNode>> new_children;
        new_children.reserve(node->children.size());
        for (auto& child : node->children) {
            if (!child) continue;
            if (same_members(node->members, child->members)) {
                for (auto& gc : child->children) {
                    if (!gc) continue;
                    gc->parent = node;
                    new_children.push_back(gc);
                }
                child->children.clear();
            } else {
                new_children.push_back(child);
            }
        }
        node->children.swap(new_children);
    }

    std::vector<int> buildFlatTreeLevel(HSTBuilder& builder,
                                        const std::vector<int>& prev_level,
                                        const HSTLevel& L) {
        std::unordered_map<int, int> center_to_parent;
        center_to_parent.reserve(L.centers.size());
        for (int c : L.centers) center_to_parent[c] = builder.addNode(c);

        const int builder_size_before = builder.size();
        std::vector<bool> attached(builder_size_before, false);
        for (int v = 0; v < n; ++v) {
            int old_node = prev_level[v];
            if (old_node >= builder_size_before || attached[old_node]) continue;
            attached[old_node] = true;
            int child_center = builder.node(old_node).center;
            int cluster = L.owner[child_center];
            builder.attach(center_to_parent.at(cluster), old_node);
        }
        for (auto& [c, idx] : center_to_parent) builder.sortMembers(idx);

        std::vector<int> new_level(n);
        for (int v = 0; v < n; ++v) {
            int child_center = builder.node(prev_level[v]).center;
            new_level[v] = center_to_parent.at(L.owner[child_center]);
        }
        return new_level;
    }

    double getMendelScalingTime() const {
        return total_time_spent_on_mendel_scaling;
    }
};

#endif //OBLIVIOUSROUTING_TREE_ORACLE_H
