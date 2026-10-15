//
// Created by Mert on 15.10.26.
//

#ifndef E_ROUTING_ROUTING_TABLE_ADAPTER_H
#define E_ROUTING_ROUTING_TABLE_ADAPTER_H
#include "allpair_routing_table.h"
#include "core/errors.h"
#include "linear_routing_table.h"
#include "routing_table_serializer.h"

#include <cstdint>


class RoutingTableAdapter {
public:
    static Result<std::uint32_t> checked_u32(int value, const char* field);

    static Result<Docu> exportLinear(const LinearRoutingTable& table, const optimized::Graph<EdgeData>& g, std::uint64_t fingerprint, std::int32_t root);
    static Result<Docu> exportAllPair(const AllPairRoutingTable& table, const optimized::Graph<EdgeData>& g, std::uint64_t fingerprint);

    static Result<void> saveLinear(const LinearRoutingTable& table, const optimized::Graph<EdgeData>& g, const std::string& path, std::uint64_t fingerprint, std::int32_t root);
    static Result<void> saveAllPair(const AllPairRoutingTable& table, const optimized::Graph<EdgeData>& g, const std::string& path, std::uint64_t fingerprint);

    static Result<LinearRoutingTable> loadLinear(const optimized::Graph<EdgeData>& g, const std::string& path, std::uint64_t fingerprint, std::int32_t* root = nullptr);
    static Result<AllPairRoutingTable> loadAllPair(const optimized::Graph<EdgeData>& g, const std::string& path, std::uint64_t fingerprint);
};


#endif // E_ROUTING_ROUTING_TABLE_ADAPTER_H
