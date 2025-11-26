//
// Created by Mert Biyikli on 24.11.25.
//

#ifndef OBLIVIOUSROUTING_EFFICIENT_CKR_H
#define OBLIVIOUSROUTING_EFFICIENT_CKR_H

#include "../../../solver/solver.h"
#include "../../../datastructures/graph_csr.h"
#include "../ckr_tree_decomposer.h"
#include "../raecke_ckr_transform.h"
#include "../utils/ultrametric_tree.h"
#include "../ckr_tree_decomposer.h"
#include <chrono>

/*
    * Optimized CKR implementation with efficient scaling:
    * - use theoretical scaling factors given by Mendel and Schwob ("Fast C-K-R Partitions of Sparse Graphs")
    * - use the precomputed MST as preprocessing
    * - embrace LCA DS for better performance and Quotient Graph Generating
 */
class EfficientCKR {

    Graph_csr m_graph;
    double m_lambdaSum = 0.0;
    double diameter = 0.0;

    MendelScaling::UltrametricTree ultrametric;

    bool use_mendel_scaling = true;      // toggle at runtime or via CLI
    uint32_t ckr_seed = 0xC0FFEE;       // for reproducible partitions


    std::vector<std::unordered_map<std::pair<int, int>, double>> edge2Load;

    std::vector<std::vector<double>> tree_edge2Load;
    public:

    // TODO:  For now these member are public however later we can make them private and provide accessors instead of exposing them
        bool debug = false;
        std::vector<double> oracle_running_times;
        std::vector<double> pure_oracle_running_times;
        // TODO: keep this for now, but can lead to memory performance issue
        std::vector<Graph_csr> m_graphs;
        std::vector<double> m_lambdas;
        std::vector<std::shared_ptr<TreeNode>> m_trees;
        void run() {
            m_lambdaSum = 0.0;
            int id = 0;
            while (m_lambdaSum < 1.0) {
                auto start = std::chrono::high_resolution_clock::now();
                m_lambdaSum += iterate(id);

                oracle_running_times.push_back(
                    std::chrono::duration<double, std::milli>(
                        std::chrono::high_resolution_clock::now() - start
                    ).count()
                );
                id++;
            }
        }


        double iterate(int id) {
            std::shared_ptr<TreeNode> t = getTree(m_graph);
            computeRLoads(t, m_graph, id);
            double l = getMaxRload(id);
            double lambda = std::min(1.0/l, 1.0 - m_lambdaSum);


            // for debuggin purposes print the tree , graph and the current edge loads, as well as the current lambda
            if (debug) {
                std::cout << "Iteration " << id << ":\n";
                std::cout << "Current lambda: " << lambda << "\n";
                std::cout << "Current Tree:\n";
                print_tree(t);
                std::cout << "Current Graph Distances:\n";
                m_graph.print();
                std::cout << "Current Edge R-Loads:\n";
                for (const auto& [edge, rLoad] : edge2Load[id]) {
                    std::cout << "Edge (" << edge.first << ", " << edge.second << ") : R-Load = " << rLoad << "\n";
                }
            }


            m_lambdas.push_back(lambda);
            m_graphs.push_back(m_graph);
            m_trees.push_back(t);
            // update weights
            m_lambdaSum += lambda;
            computeNewDistances(m_graph);

            //edge2Load.clear(); // clear for next iteration
            return lambda;
        }

        const int getIterationCount() const {
            return static_cast<int>(oracle_running_times.size());
        }


        std::vector<CKRLevel> m_levels;          // from finest (0) upward


        std::vector<int> build_ckr_level(const Graph_csr& g, double Delta, CKRLevel& L);




        void init(const Graph_csr& g);
        void preprocess();
        std::shared_ptr<TreeNode> getTree(Graph_csr& g);
        void computeRLoads(std::shared_ptr<TreeNode> t, Graph_csr& g, int tree_index);
        double getMaxRload(int tree) const;
        // void addLoadToEdge(int u, int v, double load);
        void computeNewDistances(Graph_csr& g);
        void setGraph(const Graph_csr& g);
        const Graph_csr& getGraph() const;


    };


#endif //OBLIVIOUSROUTING_EFFICIENT_CKR_H