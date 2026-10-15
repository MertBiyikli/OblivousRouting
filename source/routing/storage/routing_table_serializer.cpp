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
    if (n == 0) return false;

    if (type == RoutingTableType::Linear) {
        return id < n;
    }

    if (type == RoutingTableType::AllPair) {
        if (static_cast<uint64_t>(id) >= static_cast<uint64_t>(n) * n) return false;

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

    for (const auto& edge : docu.edges) {
        for (const auto& [id, value] : edge) {
            if (!is_valid_id(docu.type, id, docu.n)) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Invalid commodity ID: " + std::to_string(id));
            }
            if (!std::isfinite(value)) {
                return makeErrorMessage(ErrorCode::InvalidArgument, "Flow value must be finite: " + std::to_string(value));
            }
        }
    }

    return {};
}

Result<void> RoutingTableSerializer::save(const Docu& doc, const std::string& path) {
    auto res = validate(doc);
    if (!res) {
        return res;
    }

    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        return makeErrorMessage(ErrorCode::RuntimeError, "Failed to open file: " + path);
    }

    out.write("EROUTE01", 8);
    res = write_u32(out, 1);
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
    if (!in)
        throw std::runtime_error("Cannot open input file: " + path);

    char magic[8]{};
    if (!in.read(magic, 8) || std::string(magic, 8) != "EROUTE01")
        throw std::runtime_error("Invalid magic");
    if (read_u32(in) != 1)
        throw std::runtime_error("Unsupported file version");

    const auto kind_raw = read_u32(in);
    if (kind_raw != static_cast<std::uint32_t>(expected_kind))
        throw std::runtime_error("Unexpected table kind");

    Docu doc{expected_kind, read_u32(in).value(), read_u32(in).value(), read_u64(in).value(), {}};
    const std::uint64_t total = read_u64(in).value();

    if (doc.n != expected_nodes || doc.m_directed != expected_directed_edges || doc.graph_fingerprint != expected_graph_fingerprint)
        throw std::runtime_error("Graph identity mismatch");

    doc.edges.resize(doc.m_directed);
    std::uint64_t remaining = total;

    for (auto& edge : doc.edges) {
        const auto count = read_u64(in).value();
        const std::uint64_t maximum = doc.type == RoutingTableType::Linear ? std::uint64_t(doc.n) : std::uint64_t(doc.n) * (doc.n - 1) / 2;

        if (count > remaining || count > maximum)
            throw std::runtime_error("Invalid edge entry count");

        remaining -= count;
        edge.reserve(static_cast<std::size_t>(count));

        for (std::uint64_t i = 0; i < count; ++i) {
            const auto id = read_u32(in).value();
            const double value = std::bit_cast<double>(read_u64(in).value());
            edge.push_back({id, value});
        }
    }

    if (remaining != 0)
        throw std::runtime_error("Entry count mismatch");
    if (in.peek() != EOF)
        throw std::runtime_error("Trailing bytes");

    validate(doc);
    return doc;
}
