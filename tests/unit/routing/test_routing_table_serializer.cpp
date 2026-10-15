//
// Created by Mert on 15.10.26.
//


#include <catch2/catch_test_macros.hpp>
#include "routing/storage/routing_table_serializer.h"

#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>


Docu makeLinearDocument() {
    Docu doc{
        RoutingTableType::Linear,
        3,
        4,
        0x123456789ABCDEF0ULL,
        1,
        std::vector<std::vector<std::pair<std::uint32_t, double>>>(4)
    };

    doc.edges[0] = {{0, 0.5}, {2, -0.25}};
    doc.edges[1] = {{0, -0.5}, {2, 0.25}};
    return doc;
}


TEST_CASE("RoutingTableSerializer preserves fingerprint, root and signed coefficients", "[routing][serialization]") {
    const auto path = std::filesystem::temp_directory_path() / "erouting_test_v2.bin";
    auto doc = makeLinearDocument();

    auto saved = RoutingTableSerializer::save(doc, path.string());
    REQUIRE(saved.has_value());

    auto loaded = RoutingTableSerializer::load(path.string(), doc.n, doc.m_directed, doc.graph_fingerprint, RoutingTableType::Linear);
    REQUIRE(loaded.has_value());

    CHECK(loaded.value().graph_fingerprint == doc.graph_fingerprint);
    CHECK(loaded.value().root_node == 1);
    REQUIRE(loaded.value().edges.size() == doc.edges.size());

    for (std::size_t e = 0; e < doc.edges.size(); ++e) {
        REQUIRE(loaded.value().edges[e].size() == doc.edges[e].size());

        for (std::size_t i = 0; i < doc.edges[e].size(); ++i) {
            CHECK(loaded.value().edges[e][i].first == doc.edges[e][i].first);
            CHECK(std::bit_cast<std::uint64_t>(loaded.value().edges[e][i].second) == std::bit_cast<std::uint64_t>(doc.edges[e][i].second));
        }
    }

    remove(path);
}
