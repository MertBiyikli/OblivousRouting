//
// Created by Mert on 15.10.26.
//

#include "routing/storage/routing_table_serializer.h"
#include <cmath>
#include <fstream>

Result<std::uint32_t> RoutingTableSerializer::read_u32(std::istream& in) {
    u_int32_t x = 0;
    for (unsigned i = 0; i < 4; ++i) {
        const int b = in.get();
        if (b == EOF) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Truncated u32");
        }
        x |= std::uint32_t(b) << (i * 8);
    }
    return x;
}

Result<std::uint64_t> RoutingTableSerializer::read_u64(std::istream& in) {
    u_int64_t x = 0;
    for (unsigned i = 0; i < 8; ++i) {
        const int b = in.get();
        if (b == EOF) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Truncated u64");
        }
        x |= std::uint64_t(b) << (i * 8);
    }
    return x;
}

Result<void> RoutingTableSerializer::write_u32(std::ofstream& ofs, uint32_t value) {
    for (unsigned i = 0; i < 4; ++i) {
        ofs.put(static_cast<char>((value >> (i * 8)) & 255));
        if (!ofs) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Failed to write u32");
        }
    }
    return {};
}

Result<void> RoutingTableSerializer::write_u64(std::ofstream& ofs, uint64_t value) {
    for (unsigned i = 0; i < 8; ++i) {
        ofs.put(static_cast<char>((value >> (i * 8)) & 255));
        if (!ofs) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Failed to write u64");
        }
    }
    return {};
}

Result<std::uint64_t> RoutingTableSerializer::entry_count(const Docu& doc) {
    std::uint64_t total = 0;
    for (const auto& edge : doc.edges) {
        if (edge.size() > std::numeric_limits<std::uint64_t>::max() - total) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Entry count overflow");
        }
        total += edge.size();
    }
    return total;
}

bool RoutingTableSerializer::is_valid_id(const RoutingTableType& type, const uint32_t& id, const uint32_t& n) {
    if (n == 0)
        return false;

    if (type == RoutingTableType::Linear) {
        return id < n;
    }

    if (type == RoutingTableType::AllPair) {
        if (static_cast<uint64_t>(id) >= static_cast<uint64_t>(n) * n)
            return false;

        const uint32_t s = id / n;
        const uint32_t t = id % n;

        return s < t;
    }

    return false;
}

Result<void> RoutingTableSerializer::validate(const Docu& docu) {
    if (docu.n == 0) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Number of nodes cannot be zero");
    }
    if (docu.m_directed == 0) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Number of directed edges cannot be zero");
    }
    if (docu.edges.size() != docu.m_directed) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Mismatch between number of directed edges and size of edges vector");
    }

    if (docu.type == RoutingTableType::Linear) {
        if (docu.root_node < 0 || static_cast<std::uint64_t>(docu.root_node) >= docu.n) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid linear routing root");
        }
    } else if (docu.type == RoutingTableType::AllPair) {
        if (docu.root_node != -1) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "All-pairs routing table must not contain a root");
        }
    } else {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Unknown routing table type");
    }

    for (const auto& edge : docu.edges) {
        const std::uint64_t maximum = docu.type == RoutingTableType::Linear ? static_cast<std::uint64_t>(docu.n) : static_cast<std::uint64_t>(docu.n) * (docu.n - 1) / 2;

        if (edge.size() > maximum) {
            return makeErrorMessage(ErrorCode::InvalidArgument, "Too many sparse coefficients on an edge");
        }

        bool has_previous = false;
        std::uint32_t previous_id = 0;

        for (const auto& [id, value] : edge) {
            if (!is_valid_id(docu.type, id, docu.n)) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid routing coefficient ID");
            }

            if (has_previous && id <= previous_id) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Routing coefficient IDs must be strictly increasing");
            }

            if (!std::isfinite(value)) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Routing coefficient must be finite");
            }

            previous_id = id;
            has_previous = true;
        }
    }

    return {};
}

Result<void> RoutingTableSerializer::save(const Docu& doc, const std::string& path) {
    auto res = validate(doc);
    if (!res) {
        return res;
    }

    std::ofstream out(path, std::ios::binary | std::ios::out | std::ios::trunc);

    if (!out.is_open()) {
        return makeErrorMessage(ErrorCode::InvalidArgument, "Failed to open routing table for writing: " + path);
    }

    out.write("EROUTE01", 8);
    res = write_u32(out, 2);
    if (!res) {
        return res;
    }
    res = write_u32(out, static_cast<std::uint32_t>(doc.type));
    if (!res) {
        return res;
    }
    res = write_u32(out, doc.n);
    if (!res) {
        return res;
    }
    res = write_u32(out, doc.m_directed);
    if (!res) {
        return res;
    }
    res = write_u64(out, doc.graph_fingerprint);
    if (!res) {
        return res;
    }

    auto res_entry_count = entry_count(doc);
    if (!res_entry_count) {
        return getError(res_entry_count);
    }
    res = write_u64(out, res_entry_count.value());
    if (!res) {
        return res;
    }

    res = write_u32(out, std::bit_cast<std::uint32_t>(doc.root_node));
    if (!res) {
        return getError(res);
    }

    for (const auto& edge : doc.edges) {
        res = write_u64(out, edge.size());
        if (!res) {
            return res;
        }

        for (const auto& entry : edge) {
            res = write_u32(out, entry.first);
            if (!res) {
                return res;
            }
            res = write_u64(out, std::bit_cast<std::uint64_t>(entry.second));
            if (!res) {
                return res;
            }
        }
    }

    out.flush();
    if (!out) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Failed to finish writing file: " + path);
    }
    return {};
}

Result<Docu> RoutingTableSerializer::load(const std::string& path, std::uint32_t expected_nodes, std::uint32_t expected_directed_edges, std::uint64_t expected_graph_fingerprint,
                                          RoutingTableType expected_kind) {
    std::ifstream in(path, std::ios::binary);

    if (!in) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Cannot open routing table: " + path);
    }

    char magic[8]{};

    if (!in.read(magic, 8) || std::string(magic, 8) != "EROUTE01") {
        return makeErrorMessage(ErrorCode::RuntimeError, "Invalid routing table magic");
    }

    auto version = read_u32(in);
    if (!version)
        return getError(version);

    if (version.value() != 2) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Unsupported routing table version");
    }

    auto type = read_u32(in);
    if (!type)
        return getError(type);

    if (type.value() != static_cast<std::uint32_t>(expected_kind)) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Unexpected routing table type");
    }

    auto nodes = read_u32(in);
    if (!nodes)
        return getError(nodes);

    auto directed_edges = read_u32(in);
    if (!directed_edges)
        return getError(directed_edges);

    auto fingerprint = read_u64(in);
    if (!fingerprint)
        return getError(fingerprint);

    auto entry_total = read_u64(in);
    if (!entry_total)
        return getError(entry_total);

    auto root_bits = read_u32(in);
    if (!root_bits)
        return getError(root_bits);

    if (nodes.value() != expected_nodes || directed_edges.value() != expected_directed_edges) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "Routing table graph dimensions differ");
    }

    if (fingerprint.value() != expected_graph_fingerprint) {
        return makeErrorMessage(ErrorCode::InvalidGraph, "Routing table graph fingerprint mismatch");
    }

    Docu doc{expected_kind, nodes.value(), directed_edges.value(), fingerprint.value(), std::bit_cast<std::int32_t>(root_bits.value()), {}};

    const std::uint64_t maximum_per_edge = expected_kind == RoutingTableType::Linear ? static_cast<std::uint64_t>(doc.n) : static_cast<std::uint64_t>(doc.n) * (doc.n - 1) / 2;

    if (doc.n == 0 || doc.m_directed == 0 || entry_total.value() > maximum_per_edge * static_cast<std::uint64_t>(doc.m_directed)) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Invalid routing table dimensions or coefficient count");
    }

    // Check the actual file size before using header counts to allocate memory.
    const auto payload_begin = in.tellg();

    in.seekg(0, std::ios::end);
    const auto file_end = in.tellg();

    if (payload_begin < std::streampos(0) || file_end < payload_begin) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Cannot determine routing table payload size");
    }

    const auto available = static_cast<std::uint64_t>(file_end - payload_begin);

    const std::uint64_t edge_header_bytes = static_cast<std::uint64_t>(doc.m_directed) * 8;
    const std::uint64_t coefficient_bytes = entry_total.value() * 12;

    if (available != edge_header_bytes + coefficient_bytes) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Routing table payload length mismatch");
    }

    in.seekg(payload_begin);

    if (!in) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Cannot seek to routing table payload");
    }

    doc.edges.resize(doc.m_directed);

    std::uint64_t remaining = entry_total.value();

    for (auto& edge : doc.edges) {
        auto count_result = read_u64(in);
        if (!count_result)
            return getError(count_result);

        const std::uint64_t count = count_result.value();

        if (count > remaining || count > maximum_per_edge) {
            return makeErrorMessage(ErrorCode::RuntimeError, "Invalid sparse edge entry count");
        }

        remaining -= count;
        edge.reserve(static_cast<std::size_t>(count));

        bool has_previous = false;
        std::uint32_t previous_id = 0;

        for (std::uint64_t i = 0; i < count; ++i) {
            auto id_result = read_u32(in);
            if (!id_result)
                return getError(id_result);

            auto value_result = read_u64(in);
            if (!value_result)
                return getError(value_result);

            const std::uint32_t id = id_result.value();
            const double value = std::bit_cast<double>(value_result.value());

            if (!is_valid_id(doc.type, id, doc.n) || !std::isfinite(value)) {
                return makeErrorMessage(ErrorCode::RuntimeError, "Invalid routing coefficient");
            }

            if (has_previous && id <= previous_id) {
                return makeErrorMessage(ErrorCode::RuntimeError, "Unsorted or duplicate routing coefficient ID");
            }

            edge.emplace_back(id, value);
            previous_id = id;
            has_previous = true;
        }
    }

    if (remaining != 0) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Routing table entry count mismatch");
    }

    auto valid = validate(doc);
    if (!valid)
        return getError(valid);

    return doc;
}

Result<RoutingTableHeader> RoutingTableSerializer::inspect(const std::string& path) {
    std::ifstream in(path, std::ios::binary);

    if (!in) {
        return makeErrorMessage(ErrorCode::FileNotFound, "Cannot open routing table: " + path);
    }

    char magic[8]{};

    if (!in.read(magic, 8) || std::string(magic, 8) != "EROUTE01") {
        return makeErrorMessage(ErrorCode::RuntimeError, "Invalid routing-table magic");
    }

    auto version = read_u32(in);
    if (!version)
        return getError(version);

    if (version.value() != 2) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Unsupported routing-table version");
    }

    auto kind = read_u32(in);
    if (!kind)
        return getError(kind);

    if (kind.value() > static_cast<std::uint32_t>(RoutingTableType::AllPair)) {
        return makeErrorMessage(ErrorCode::InvalidRouting, "Unknown routing-table representation");
    }

    return RoutingTableHeader{static_cast<RoutingTableType>(kind.value()), version.value()};
}
