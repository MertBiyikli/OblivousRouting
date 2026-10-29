//
// Created by Mert on 29.10.26.
//

#include "analysis/traffic_change/demand_vector.h"

#include <cmath>
#include <limits>
#include <map>
#include <utility>

Result<DemandVector> DemandVector::create(
    int nodeCount,
    const std::vector<DemandEntry>& entries
) {
    if (nodeCount <= 0) {
        return makeErrorMessage(
            ErrorCode::InvalidDemand,
            "Node count must be positive"
        );
    }

    std::map<std::pair<int, int>, double> aggregated;

    for (const auto& entry : entries) {
        if (entry.source < 0 || entry.source >= nodeCount ||
            entry.target < 0 || entry.target >= nodeCount ||
            entry.source == entry.target) {
            return makeErrorMessage(
                ErrorCode::InvalidDemand,
                "Invalid demand endpoints"
            );
        }

        if (!std::isfinite(entry.volume) || entry.volume < 0.0) {
            return makeErrorMessage(
                ErrorCode::InvalidDemand,
                "Demand volume must be finite and nonnegative"
            );
        }

        // Generated demand models can produce zero-demand pairs.
        if (entry.volume == 0.0) {
            continue;
        }

        auto& sum = aggregated[{entry.source, entry.target}];
        sum += entry.volume;

        if (!std::isfinite(sum)) {
            return makeErrorMessage(
                ErrorCode::InvalidDemand,
                "Commodity demand overflow"
            );
        }
    }

    if (aggregated.empty()) {
        return makeErrorMessage(
            ErrorCode::InvalidDemand,
            "Demand vector contains no positive traffic"
        );
    }

    DemandVector result;
    result.node_count_ = nodeCount;
    result.entries_.reserve(aggregated.size());

    for (const auto& [pair, volume] : aggregated) {
        result.entries_.push_back({
            pair.first,
            pair.second,
            volume
        });

        result.total_volume_ += volume;

        if (!std::isfinite(result.total_volume_)) {
            return makeErrorMessage(
                ErrorCode::InvalidDemand,
                "Total demand volume overflow"
            );
        }
    }

    return result;
}

Result<DemandVector> DemandVector::fromDemands(
    const demands& input,
    int nodeCount
) {
    // Check lengths directly; demands::size() asserts on mismatch.
    const auto size = input.demand_values.size();

    if (input.source.size() != size ||
        input.target.size() != size) {
        return makeErrorMessage(
            ErrorCode::InvalidDemand,
            "Inconsistent engine demand arrays"
        );
    }

    std::vector<DemandEntry> entries;
    entries.reserve(size);

    for (std::size_t i = 0; i < size; ++i) {
        entries.push_back({
            input.source[i],
            input.target[i],
            input.demand_values[i]
        });
    }

    return create(nodeCount, entries);
}

Result<DemandVector> DemandVector::scaled(double factor) const {
    if (!std::isfinite(factor) || factor <= 0.0) {
        return makeErrorMessage(
            ErrorCode::InvalidDemand,
            "Scale factor must be finite and positive"
        );
    }

    std::vector<DemandEntry> transformed = entries_;

    for (auto& entry : transformed) {
        entry.volume *= factor;

        if (!std::isfinite(entry.volume) ||
            entry.volume <= 0.0) {
            return makeErrorMessage(
                ErrorCode::InvalidDemand,
                "Scaled demand overflow or underflow"
            );
        }
    }

    return create(node_count_, transformed);
}

Result<DemandVector> DemandVector::scaledCommodity(
    int source,
    int target,
    double factor
) const {
    if (source < 0 || source >= node_count_ ||
        target < 0 || target >= node_count_ ||
        source == target ||
        !std::isfinite(factor) || factor <= 0.0) {
        return makeErrorMessage(
            ErrorCode::InvalidDemand,
            "Invalid commodity transformation"
        );
    }

    std::vector<DemandEntry> transformed = entries_;
    bool found = false;

    for (auto& entry : transformed) {
        if (entry.source == source &&
            entry.target == target) {
            entry.volume *= factor;
            found = true;

            if (!std::isfinite(entry.volume) ||
                entry.volume <= 0.0) {
                return makeErrorMessage(
                    ErrorCode::InvalidDemand,
                    "Commodity demand overflow or underflow"
                );
            }

            break;
        }
    }

    if (!found) {
        return makeErrorMessage(
            ErrorCode::InvalidDemand,
            "Commodity does not exist in demand vector"
        );
    }

    return create(node_count_, transformed);
}

Result<demands> DemandVector::toDemands() const {
    demands output;

    for (const auto& entry : entries_) {
        // Every normalized directed pair appears exactly once.
        output.addDemand(
            entry.source,
            entry.target,
            entry.volume
        );
    }

    return output;
}

const std::vector<DemandEntry>&
DemandVector::entries() const noexcept {
    return entries_;
}

double DemandVector::totalVolume() const noexcept {
    return total_volume_;
}

int DemandVector::nodeCount() const noexcept {
    return node_count_;
}
