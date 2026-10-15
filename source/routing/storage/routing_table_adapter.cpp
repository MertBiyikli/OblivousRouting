//
// Created by Mert on 15.10.26.
//

#include <bitset>
#include "routing/storage/routing_table_adapter.h"
#include "routing/storage/routing_table_serializer.h"

Result<uint32_t> RoutingTableAdapter::checked_u32(int value, const char* field) {
    if (value < 0 || value > std::numeric_limits<std::uint32_t>::max()) {
        return makeErrorMessage(ErrorCode::InvalidArgument, std::string(field) + " is out of range for u32");
    }
    return static_cast<std::uint32_t>(value);
}

Result<Docu> RoutingTableAdapter::exportLinear(const LinearRoutingTable& table, const optimized::Graph<EdgeData>& g, u_int32_t fingerprint) {
    auto m = checked_u32(g.getNumDirectedEdges(), "edge_count(directed)");
    if (!m) {
        return getError(m);
    }

    auto n = checked_u32(g.getNumNodes(), "node_count");
    if (!n) {
        return getError(n);
    }

    if (table.n != g.getNumNodes() || table.src_ids.size() != m || table.src_flows.size() != m) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Linear table dimensions differ from graph");
    }


    Docu doc{RoutingTableType::Linear,
        n.value(),
        m.value(),
        fingerprint,
        table.getRootX(),
        std::vector<std::vector<std::pair<u_int32_t, double>>>(m.value())};

    for (std::size_t e = 0; e < m.value(); ++e) {
        if (table.src_ids[e].size() != table.src_flows[e].size()) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Mismatched linear coefficient arrays");
        }

        for (std::size_t i = 0; i < table.src_ids[e].size(); ++i) {
            auto check = checked_u32(table.src_ids[e][i], "source ID");
            if (!check) {
                return getError(check);
            }
            doc.edges[e].push_back({check.value(), table.src_flows[e][i]});
        }
    }

    auto result = RoutingTableSerializer::validate(doc);
    if (!result) {
        return getError(result);
    }
    return doc;
}

Result<Docu> RoutingTableAdapter::exportAllPair(const AllPairRoutingTable& table, const optimized::Graph<EdgeData>& g, u_int64_t fingerprint) {
    auto m = checked_u32(g.getNumDirectedEdges(), "edge_count(directed)");
    if (!m) {
        return getError(m);
    }

    auto n = checked_u32(g.getNumNodes(), "node_count");
    if (!n) {
        return getError(n);
    }

    if (table.n != g.getNumNodes() || table.adj_ids.size() != m || table.adj_vals.size() != m) {
        throw std::runtime_error("All-pair table dimensions differ from graph");
    }

    Docu doc{RoutingTableType::AllPair,
        n.value(),
        m.value(),
        fingerprint,
        -1,
        std::vector<std::vector<std::pair<u_int32_t, double>>>(m.value())};

    for (std::size_t e = 0; e < m.value(); ++e) {
        if (table.adj_ids[e].size() != table.adj_vals[e].size()) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Mismatched all-pair coefficient arrays");
        }

        for (std::size_t i = 0; i < table.adj_ids[e].size(); ++i) {
            auto check = checked_u32(table.adj_ids[e][i], "commodity ID");
            if (!check) {
                return getError(check);
            }
            doc.edges[e].push_back({check.value(), table.adj_vals[e][i]});
        }
    }

    auto result = RoutingTableSerializer::validate(doc);
    if (!result) {
        return getError(result);
    }
    return doc;
}

Result<void> RoutingTableAdapter::saveLinear(const LinearRoutingTable& table, const optimized::Graph<EdgeData>& g, const std::string& path, const u_int32_t& fingerprint) {
    auto docu = exportLinear(table, g, fingerprint);
    if (!docu) {
        return getError(docu);
    }
    return RoutingTableSerializer::save(docu.value(), path);
}

Result<void> RoutingTableAdapter::saveAllPair(const AllPairRoutingTable& table, const optimized::Graph<EdgeData>& g, const std::string& path, const u_int64_t& fingerprint) {
    auto docu = exportAllPair(table, g, fingerprint);
    if (!docu) {
        return getError(docu);
    }
    return RoutingTableSerializer::save(docu.value(), path);
}

Result<LinearRoutingTable> RoutingTableAdapter::loadLinear(const optimized::Graph<EdgeData>& g, const std::string& path, const u_int32_t& fingerprint) {
    auto check_n = checked_u32(g.getNumNodes(), "node count");
    if (!check_n) {
        return getError(check_n);
    }
    auto check_directed_edges = checked_u32(g.getNumDirectedEdges(), "directed edge count");
    if (!check_directed_edges) {
        return getError(check_directed_edges);
    }

    auto doc = RoutingTableSerializer::load(path,
        check_n.value(),
        check_directed_edges.value(),
        fingerprint, RoutingTableType::Linear);

    if (!doc) {
        return getError(doc);
    }

    LinearRoutingTable table;
    table.init(g);

    for (std::size_t e = 0; e < doc.value().edges.size(); ++e) {
        for (const auto& entry : doc.value().edges[e]) {
            if (entry.first > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
                return makeErrorMessage(ErrorCode::RuntimeError, "Source ID exceeds int range");
            }
            table.src_ids[e].push_back(static_cast<int>(entry.first));
            table.src_flows[e].push_back(entry.second);
        }
    }

    return table;
}




Result<AllPairRoutingTable> RoutingTableAdapter::loadAllPair(const optimized::Graph<EdgeData>& g, const std::string& path, const u_int64_t& fingerprint) {
    auto check_n = checked_u32(g.getNumNodes(), "node count");
    if (!check_n) {
        return getError(check_n);
    }
    auto check_directed_edges = checked_u32(g.getNumDirectedEdges(), "directed edge count");
    if (!check_directed_edges) {
        return getError(check_directed_edges);
    }

    auto doc = RoutingTableSerializer::load(path,
        check_n.value(),
        check_directed_edges.value(),
        fingerprint, RoutingTableType::AllPair);

    if (!doc) {
        return getError(doc);
    }

    AllPairRoutingTable table;
    table.init(g);

    for (std::size_t e = 0; e < doc.value().edges.size(); ++e) {
        for (const auto& entry : doc.value().edges[e]) {
            if (entry.first > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
                return makeErrorMessage(ErrorCode::RuntimeError, "Commodity ID exceeds int range");
            }
            table.adj_ids[e].push_back(static_cast<int>(entry.first));
            table.adj_vals[e].push_back(entry.second);
        }
    }

    return table;
}
