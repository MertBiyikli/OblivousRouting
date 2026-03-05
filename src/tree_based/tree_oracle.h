//
// Created by Mert Biyikli on 09.02.26.
//

#ifndef OBLIVIOUSROUTING_TREE_ORACLE_H
#define OBLIVIOUSROUTING_TREE_ORACLE_H

#include <list>

#include "hst.h"
#include "../datastructures/IGraph.h"
#include "utils/ultrametric_tree.h"
#include "utils/quotient_graph.h"
#include "mst.h"
#include <vector>

class TreeOracle {
public:
    explicit TreeOracle(IGraph& graph) : graph(graph) {
        n=graph.getNumNodes();
        diameter=graph.GetDiameter();
    }

    TreeOracle(IGraph& graph, bool activateMendelScaling) : graph(graph) {
        n=graph.getNumNodes();
        diameter=graph.GetDiameter();
        this->applyMendelScaling = activateMendelScaling;
        total_time_spent_on_mendel_scaling = 0.0;
    }
    virtual ~TreeOracle() = default;
    IGraph& graph;
    int n;
    double diameter;
    std::vector<double> scales; // the scales at which to compute the partitions
    std::vector<int> perm;

    bool applyMendelScaling = false; // whether to apply Mendel scaling to the distances before computing the tree
    MendelScaling::UltrametricTree ultrametric;
    MendelScaling::QuotientConstruction qc;
    double total_time_spent_on_mendel_scaling;

    virtual std::shared_ptr<HSTNode> getTree(std::vector<double>& distances) {
        total_time_spent_on_mendel_scaling = 0;
        updateDistances(distances);
        std::shared_ptr<HSTNode> tree = computeHST();
        return tree;
    }

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

        std::shared_ptr<HSTNode> root =  std::make_shared<HSTNode>(std::numeric_limits<int>::max());

        if (applyMendelScaling) {
            auto start = timeNow();
            qc.preprocessEdges(graph);
            total_time_spent_on_mendel_scaling += duration(timeNow() - start);
        }


        std::vector<std::shared_ptr<HSTNode>> current_tree_nodes = prev_nodes; // start with the leaf nodes


        for (const double Delta : scales) {
            std::vector<std::shared_ptr<HSTNode>> current_level;
            HSTLevel L;


            MendelScaling::QuotientLevel Q;
            if (applyMendelScaling) {
                auto start = timeNow();
                Q = qc.constructQuotientGraph(ultrametric, Delta, graph);
                //Q = build_quotient_graph_using_sliding_window(graph, ultrametric, Delta);
                if (Q.Gq->getNumNodes() <= 1) {
                    // if single cluster: make that the root once at the end
                    continue;
                }

                total_time_spent_on_mendel_scaling += duration(timeNow() - start);
                computeQuotientLevelPartition(Q, L, Delta);

            }else {
                // compute the partition on the original graph and build the tree level
                computeLevelPartition(graph, L, perm, Delta);
            }
            current_level = buildTreeLevel(prev_nodes, L);


            for (auto& node : current_level) {
                current_tree_nodes.push_back(node);
            }

            prev_nodes = current_level;
        }
        finishTree(root, current_tree_nodes);
        return root;
    }

    void preprocess() {
        n = graph.getNumNodes();
        if (n == 0) throw std::invalid_argument("The graph has no nodes.");
        diameter = graph.GetDiameter();

        perm.clear();
        // initialize the node permutation
        for (int v = 0; v < n; ++v) {
            perm.push_back(v);
        }

        if (applyMendelScaling) {
            auto start = timeNow();
            RandomMST mst_oracle(graph);
            auto mst = mst_oracle.computeMST();

            std::vector<double> mst_weights;
            mst_weights.reserve(mst.size());
            for (auto [u,v] : mst) {
                mst_weights.push_back(graph.getEdgeDistance(u,v));
            }

            ultrametric.buildFromMST(graph.getNumNodes(), mst, mst_weights);
            total_time_spent_on_mendel_scaling += duration(timeNow() - start);
            assert(ultrametric.root != -1);
        }
    }


    // Helper: build one representative original vertex per quotient node.
    // sigma[v] = quotient id of original vertex v  (size nG, values in [0..nQ-1])
    static std::vector<int> buildRepresentatives(const std::vector<int>& sigma, int nQ) {
        const int nG = (int)sigma.size();
        std::vector<int> rep(nQ, -1);
        for (int v = 0; v < nG; ++v) {
            int q = sigma[v];
            if (q < 0 || q >= nQ) {
                throw std::runtime_error("sigma[v] out of range at v=" + std::to_string(v));
            }
            if (rep[q] == -1) rep[q] = v;
        }
        for (int q = 0; q < nQ; ++q) {
            if (rep[q] == -1) {
                throw std::runtime_error("Quotient node q=" + std::to_string(q) + " has no preimage (rep=-1)");
            }
        }
        return rep;
    }

    void computeQuotientLevelPartition(MendelScaling::QuotientLevel& Q, HSTLevel& level, double delta) {
        // 1) Compute partition on quotient graph
        HSTLevel qL;
        std::vector<int> x_perm;
        std::vector<char> used(Q.Gq->getNumNodes(), 0);

        for (auto v : perm) {
            int vq = Q.sigma_compact_of_v[v];
            if (vq < 0 || vq >= Q.Gq->getNumNodes()) {
                throw std::runtime_error("computeQuotientLevelPartition: sigma_compact_of_v[v] out of range at v=" + std::to_string(v));
            }
            if (!used[vq]) {
                x_perm.push_back(vq);
                used[vq] = 1;
            }
        }

        computeLevelPartition(*Q.Gq, qL, x_perm, delta);

        const int nG = graph.getNumNodes();

        // 2) Map quotient partition back to original graph
        level.R = qL.R;
        level.centers.clear();
        level.owner.resize(nG);
        for (int v = 0; v < nG; ++v) {
                level.owner[v] =v; // initialize owner to -1 (unassigned)
            }

        // identify the cluster centers for Q
        std::unordered_map<int, int> q_center_to_g_center;
        for (size_t c = 0; c < qL.centers.size(); ++c) {
            int q_center = qL.centers[c];
            int g_center = Q.members_of_q[q_center].empty() ? throw std::runtime_error("computeQuotientLevelPartition: node cluster empty!") : Q.members_of_q[q_center][0];
            level.centers.push_back(g_center);
            q_center_to_g_center[q_center] = g_center;
        }
        for (int vq = 0; vq<Q.Gq->getNumNodes(); ++vq)  {
            int cluster_q = qL.owner[vq];
            if (cluster_q == -1) {
                throw std::runtime_error("computeQuotientLevelPartition: node " + std::to_string(vq) + " in quotient graph has no cluster!");
            }
            int cluster_g = q_center_to_g_center[cluster_q];
            for (int v : Q.members_of_q[vq]) {
                level.owner[v] = cluster_g;
            }
        }
    }


    virtual void computeLevelPartition(IGraph& g, HSTLevel& level, const std::vector<int>& x_perm, double delta) = 0;


    void updateDistances(std::vector<double>& distances) {
        for (int e = 0; e < graph.getNumEdges(); ++e) {
            graph.updateEdgeDistance(e, distances[e]);
        }
    }

    void computeScales() {
        scales.clear();
        if (n == 0) throw std::invalid_argument("TreeOracle: The graph has no nodes.");
        if (applyMendelScaling){
            computeMendelScales();
        } else {
            computeNaiveScales();
        }
    }

    std::vector<std::shared_ptr<HSTNode>> buildTreeLevel(
        std::vector<std::shared_ptr<HSTNode>>& prev_nodes,
        const HSTLevel& L) {

        // Build tree level based on the partition L
        std::vector<std::shared_ptr<HSTNode>> parents(L.centers.size());

        for (size_t c = 0; c < L.centers.size(); ++c) {
            const int center = L.centers[c];
            parents[c] = std::make_shared<HSTNode>(-1);
            parents[c]->center = center;
            parents[c]->parent.reset();
            parents[c]->children.clear();

        }

        // Attach lower-level nodes to parents according to the partition
        for (const auto& child : prev_nodes) {
            if (!child) continue;

            int child_center = child->center;
            int cluster = L.owner[child_center];
            // find the according parent
            int c_id = -1;
            for (size_t i = 0; i < L.centers.size(); ++i) {
                if (L.centers[i] == cluster) {
                    c_id = i;
                    break;
                }
            }
            auto& parent = parents[c_id];

            // make them a family :)
            child->parent = parent;
            parent->children.push_back(child);

            parent->members.insert(
                        parent->members.end(),
                        child->members.begin(),
                        child->members.end()
                    );
        }

        for (auto& p : parents) {
            auto& M = p->members;
            std::sort(M.begin(), M.end());
            M.erase(std::unique(M.begin(), M.end()), M.end());
        }

        // Clean empty nodes
        std::vector<std::shared_ptr<HSTNode>> cleaned_parents;
        for (int i = 0; i < parents.size(); ++i) {
            auto& p = parents[i];
            if (p->members.size() != 0) {
                cleaned_parents.push_back(p);
            }
        }

        return cleaned_parents;
    }

    void computeMendelScales() {
        double min_weight = std::numeric_limits<double>::max();
        for (int e = 0; e< graph.getNumEdges(); e++) {
            double w = graph.getEdgeDistance(e);
            if (w < min_weight) {
                min_weight = w;
            }
        }

        // compute the new diameter after preprocessing
        if (diameter == 0) {
            throw std::runtime_error("TreeOracle: Graph has zero diameter; cannot build CKR tree.");
        }
        double d = diameter * 2 * n;
        for (; d >= diameter / (2.0 * n) || d > min_weight; d /= 8.0)
            scales.push_back(d);

        // sort the scales in increasing order
        std::reverse(scales.begin(), scales.end());
    }

    void computeNaiveScales() {
        int i = static_cast<int>(std::ceil(std::log2(diameter) / std::log2(2.0))) + 1;
        for (; i >= 0; --i) {
            double scale = (double)(1ull << i);
            if (scale <= diameter) {
                scales.push_back(scale);
            }
        }

        std::reverse(scales.begin(), scales.end());
    }

    void finishTree(std::shared_ptr<HSTNode>& root, std::vector<std::shared_ptr<HSTNode>>& prev_nodes) {
        if (prev_nodes.empty()) {
            root = std::make_shared<HSTNode>(std::numeric_limits<int>::max());
            root->members.resize(graph.getNumNodes());
            std::iota(root->members.begin(), root->members.end(), 0);
            root->center = (!root->members.empty() ? root->members[0] : -1 );
        } else {
            root = std::make_shared<HSTNode>(std::numeric_limits<int>::max());
            root->members.resize(graph.getNumNodes());
            std::iota(root->members.begin(), root->members.end(), 0);
            root->center = prev_nodes.back()->center;

            std::set<int> members_set;
            for (int i = prev_nodes.size()-1; i >= 0; --i) {
                auto& node = prev_nodes[i];
                for (auto& v : node->members) {
                    members_set.insert(v);
                }
                root->children.push_back(node);
                if (members_set.size() == graph.getNumNodes()) {
                    break;
                }
            }
        }

        cleanUpTree(root);
    }
    void cleanUpTree(std::shared_ptr<HSTNode>& node) {
        if (!node) return;

        if (node->members.size() == 0) {
            // drop empty nodes
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

        // first clean children
        for (auto& ch : node->children) cleanUpTree(ch);

        // then rebuild children list safely
        std::vector<std::shared_ptr<HSTNode>> new_children;
        new_children.reserve(node->children.size());

        for (auto& child : node->children) {
            if (!child) continue;

            if (same_members(node->members, child->members)) {
                // absorb child: adopt grandchildren
                for (auto& gc : child->children) {
                    if (!gc) continue;
                    gc->parent = node;
                    new_children.push_back(gc);
                }
                child->children.clear();
                // child gets dropped
            } else {
                new_children.push_back(child);
            }
        }

        node->children.swap(new_children);
    }

    double getMendelScalingTime() const {
        return total_time_spent_on_mendel_scaling;
    }
};

#endif //OBLIVIOUSROUTING_TREE_ORACLE_H