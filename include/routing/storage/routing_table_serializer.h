//
// Created by Mert on 15.10.26.
//

#ifndef E_ROUTING_ROUTING_TABLE_SERIALIZER_H
#define E_ROUTING_ROUTING_TABLE_SERIALIZER_H

#include "routing/routing_table.h"
#include <core/errors.h>
#include <string>

enum class RoutingTableType {
    Linear = 0,
    AllPair = 1
};

struct Docu {
    RoutingTableType type;
    u_int32_t n;
    u_int32_t m_directed;
    u_int64_t graph_fingerprint;
    std::vector<std::vector<std::pair<u_int32_t, double>>> edges;
};

class RoutingTableSerializer {
    // helper functions for serialization and deserialization of routing tables
    static Result<std::uint32_t> read_u32(std::istream& in);
    static Result<std::uint64_t> read_u64(std::istream& in);

    static Result<void> write_u32(std::ofstream& ofs, uint32_t value);
    static Result<void> write_u64(std::ofstream& ofs, uint64_t value);

    static bool is_valid_id(const RoutingTableType& type, const u_int32_t& id, const u_int32_t& n);
    static Result<std::uint64_t> entry_count(const Docu& doc);
public:
    static Result<void> validate(const Docu& docu);
    static Result<void> save(const Docu& table, const std::string& path);
    static Result<Docu> load(const std::string& path, std::uint32_t expected_nodes, std::uint32_t expected_directed_edges,
                     std::uint64_t expected_graph_fingerprint, RoutingTableType expected_kind);


};


#endif // E_ROUTING_ROUTING_TABLE_SERIALIZER_H
