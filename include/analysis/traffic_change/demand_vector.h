//
// Created by Mert on 29.10.26.
//

#ifndef E_ROUTING_DEMAND_VECTOR_H
#define E_ROUTING_DEMAND_VECTOR_H
#include "core/errors.h"
#include "utils/demands.h"

struct DemandEntry {
    int source;
    int target;
    double volume;
};

class DemandVector {
  public:
    static Result<DemandVector> fromDemands(const demands& input, int nodeCount);

    static Result<DemandVector> create(int nodeCount, const std::vector<DemandEntry>& entries);

    Result<DemandVector> scaled(double factor) const;

    Result<DemandVector> scaledCommodity(int source, int target, double factor) const;

    Result<demands> toDemands() const;

    const std::vector<DemandEntry>& entries() const noexcept;
    double totalVolume() const noexcept;
    int nodeCount() const noexcept;

  private:
    DemandVector() = default;

    int node_count_ = 0;
    double total_volume_ = 0.0;
    std::vector<DemandEntry> entries_;
};

#endif // E_ROUTING_DEMAND_VECTOR_H
