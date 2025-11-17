//
// Created by Mert Biyikli on 17.11.25.
//

#include <gtest/gtest.h>
#include <chrono>
#include "../../src/graph_csr.h"
#include "../../src/graph.h"


template<typename T>
std::vector<double> measureEdgeAdditionTime(T& graph, int const n,const int m) {
    std::vector<double> run_times;
    run_times.reserve(m);

    std::chrono::high_resolution_clock::time_point begin, end;

    int ctr_added_edges = 0;
    for (int i = 0; i<n; ++i) {
        for (int j = i+1; j<n; ++j) {

            if ( ctr_added_edges >= m ) break;
            try {
                begin = std::chrono::high_resolution_clock::now();
                graph.addEdge(i, j, 1.0);
                end = std::chrono::high_resolution_clock::now();
                run_times.push_back(std::chrono::duration<double, std::milli>(end - begin).count());
                ctr_added_edges++;
            } catch (const std::out_of_range& e) {
                // Ignore out of range errors
            }
        }
    }


    return run_times;
}


// First test for testing the scalability of adding edges to the graph. Measure running time for adding 1 million edges.
TEST(StressTest, CSR_LargeGraphEdgeAddition) {
    constexpr int num_nodes = 100000;
    constexpr int num_edges = 1000000;
    Graph_csr g(num_nodes);

    std::vector<double> run_times = measureEdgeAdditionTime(g, num_nodes, num_edges);
    double total_time = 0.0;
    for (const auto& t : run_times) {
        total_time += t;
    }
    std::cout << "Average time for adding " << num_edges << " random nodes for edges: " << total_time/run_times.size() << " ms\n";
}


// First test for testing the scalability of adding edges to the graph. Measure running time for adding 1 million edges.
TEST(StressTest, LargeGraphEdgeAddition) {
    const int num_nodes = 100000;
    const int num_edges = 1000000;
    Graph g(num_nodes);

    std::vector<double> run_times = measureEdgeAdditionTime(g, num_nodes, num_edges);
    double total_time = 0.0;
    for (const auto& t : run_times) {
        total_time += t;
    }
    std::cout << "Average time for adding " << num_edges << " random nodes for edges: " << total_time/run_times.size() << " ms\n";
}



/* --------- Edge Weight Update Stress Tests --------- */
TEST(StressTest, CSR_LargeGraphEdgeUpdate) {
    constexpr int num_nodes = 10000;
    constexpr int num_edges = 100000;
    constexpr int num_updates = 1000000;
    Graph_csr g(num_nodes);

    // First add edges
    measureEdgeAdditionTime(g, num_nodes, num_edges);

    g.preprocess();


    // Now measure update times
    std::vector<double> run_times;
    run_times.reserve(num_edges);

    std::chrono::high_resolution_clock::time_point begin, end;

    int ctr_updated_edges = 0;

    while (ctr_updated_edges < num_updates) {
        for (int i = 0; i<num_nodes; ++i) {
            for (int j = i+1; j<num_nodes; ++j) {

                if ( ctr_updated_edges >= num_updates ) break;
                try {
                    begin = std::chrono::high_resolution_clock::now();
                    g.updateEdgeDistance(i, j, 2.0);
                    end = std::chrono::high_resolution_clock::now();
                    run_times.push_back(std::chrono::duration<double, std::milli>(end - begin).count());
                    ctr_updated_edges++;
                } catch (const std::runtime_error& e) {
                    // Ignore out of range errors
                }
            }
        }
    }


    double total_time = 0.0;
    for (const auto& t : run_times) {
        total_time += t;
    }
    std::cout << "Average time for updating " << num_edges << " random edges: " << total_time/run_times.size() << " ms\n";
}


TEST(StressTest, LargeGraphEdgeUpdate) {
    const int num_nodes = 10000;
    const int num_edges = 100000;
    constexpr int num_updates = 1000000;
    Graph g(num_nodes);

    // First add edges
    measureEdgeAdditionTime(g, num_nodes, num_edges);

    // Now measure update times
    std::vector<double> run_times;
    run_times.reserve(num_edges);

    std::chrono::high_resolution_clock::time_point begin, end;

    int ctr_updated_edges = 0;
    while (ctr_updated_edges < num_updates) {
        for (int i = 0; i<num_nodes; ++i) {
            for (int j = i+1; j<num_nodes; ++j) {

                if ( ctr_updated_edges >= num_updates ) break;
                try {
                    begin = std::chrono::high_resolution_clock::now();
                    g.updateEdgeDistance(i, j, 2.0);
                    end = std::chrono::high_resolution_clock::now();
                    run_times.push_back(std::chrono::duration<double, std::milli>(end - begin).count());
                    ctr_updated_edges++;
                } catch (const std::runtime_error& e) {
                    // Ignore out of range errors
                }
            }
        }
    }

    double total_time = 0.0;
    for (const auto& t : run_times) {
        total_time += t;
    }
    std::cout << "Average time for updating " << num_edges << " random edges: " << total_time/run_times.size() << " ms\n";
}