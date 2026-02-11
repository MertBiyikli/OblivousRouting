//
// Created by Mert Biyikli on 07.02.26.
//

#include "xcut_expander_hierarchy.h"

#include <iostream>

#include "../datastructures/IGraph.h"
#include <utility>
#include <vector>

#include "xcut/expanders/expander_hierarchy.hpp"
#include "xcut/expanders/expander_decomposition.hpp"
#include "xcut/expanders/refinement.hpp"
#include "xcut/data_structures/graph.hpp"
#include "xcut/core/definitions.hpp"
#include "xcut/core/config.hpp"
#include "xcut/algorithms/bfs.hpp"


#include "spdlog/async.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/ranges.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/spdlog.h"



inline void ensure_xcut_logger(bool verbose=false, bool debug=false) {
    if (spdlog::get("xcut")) return;

    auto null_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    auto logger = std::make_shared<spdlog::logger>("xcut", null_sink);

    if (debug)       logger->set_level(spdlog::level::debug);
    else if (verbose)logger->set_level(spdlog::level::info);
    else             logger->set_level(spdlog::level::warn);

    spdlog::register_logger(logger);
}


void ExpanderRouting::init(IGraph &g) {
    const bool verbose = true;
    const bool debug   = true;
    ensure_xcut_logger(verbose, debug);
    auto xcut_logger = spdlog::get("xcut");
    std::vector<Edge> edges;
    std::vector<EdgeWeight> weights;
    edges.reserve(g.getNumEdges());
    weights.reserve(g.getNumEdges());

    constexpr double SCALE = 100000.0; // pick 1e4..1e6 depending on your cap range

    for (int e = 0; e < g.getNumEdges(); ++e) {
        auto [u, v] = g.edgeEndpoints(e);

        // skip self-loops if any
        if (u == v) continue;

        // keep each undirected edge once
        if (u > v) continue;

        // If your IGraph stores both directions as separate edges,
        // you must avoid duplicates here. Easiest: only accept when u < v AND
        // additionally only if (e corresponds to the "canonical" direction).
        // If you don't have that, use a hash set (see note below).
        // For now, accept u<v but you'll still get duplicates if both arcs exist.
        // --> Better: use a set, see the note.

        double cap = g.getEdgeCapacity(e);
        unsigned w = (unsigned) llround(cap * SCALE);
        if (w == 0) w = 1;

        edges.push_back(Edge{(NodeID)u, (NodeID)v});
        weights.push_back((EdgeWeight)w);
    }

    int num_clusters = 2; // Example value, set as needed
    auto config = Config(num_clusters);
    config.m_verbose = true;
    config.m_debug = true;


    // 2) Build xcut graph as UNDIRECTED (directed=false)
    Graph xg(edges, weights, /*directed=*/false);

    if (xg.has_degree_zero()) {
        std::cerr << "Error: Input graph has degree-zero vertices.\n";
        return;
    }
    // 4) Run the exact pipeline from xcut main.cpp
    auto sparsifier = expander_hierarchy(&xg, &config);

    (void) normalized_cut(sparsifier, config.m_num_clusters);

    for (NodeID level = sparsifier.size(); level > 0; --level) {
        refine(sparsifier, level - 1, config);
    }


    auto clustering = sparsifier.clustering();
    xcut_logger->info("clustering: {}", fmt::join(clustering, ", "));

    std::vector<int> cnt(config.m_num_clusters, 0);
    for (auto c : clustering) cnt[c]++;
    xcut_logger->info("cluster sizes: {}", fmt::join(cnt, ", "));


}
