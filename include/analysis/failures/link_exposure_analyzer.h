//
// Created by Mert on 29.10.26.
//

#ifndef E_ROUTING_LINK_EXPOSURE_ANALYZER_H
#define E_ROUTING_LINK_EXPOSURE_ANALYZER_H
#include "core/errors.h"
#include "routing/routing_table.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

struct LinkExposureData {
    int edge_id = -1;
    int source = -1;
    int target = -1;

    double affected_demand_volume = 0.0;
    double lost_traffic = 0.0;
    double lost_traffic_fraction = 0.0;
};

struct LinkExposureSummary {
    std::size_t tested_links = 0;

    int most_critical_edge_id = -1;
    int most_critical_source = -1;
    int most_critical_target = -1;

    double maximum_lost_traffic_fraction = 0.0;
    double average_lost_traffic_fraction = 0.0;
    double median_lost_traffic_fraction = 0.0;

    double maximum_affected_demand_fraction = 0.0;

    std::size_t traffic_carrying_links = 0;

    std::size_t critical_links_10_percent = 0;
    std::size_t critical_links_25_percent = 0;
    std::size_t critical_links_50_percent = 0;

    std::vector<LinkExposureData> links;
};

class LinkExposureAnalyzer {
  public:
    static Result<LinkExposureSummary> analyze(const optimized::Graph<EdgeData>& graph, const RoutingScheme& scheme, const demands& demand) {
        constexpr double flow_epsilon = 1e-12;

        double total_demand = 0.0;

        for (std::size_t i = 0; i < demand.size(); ++i) {
            const auto [s, t] = demand.getDemandPair(i);
            const double value = demand.getDemandValue(i);

            if (s < 0 || t < 0 || s >= graph.getNumNodes() || t >= graph.getNumNodes() || !std::isfinite(value) || value < 0.0) {
                return makeErrorMessage(ErrorCode::InvalidDemand, "Invalid commodity in link exposure analysis.");
            }

            if (s == t && value > 0.0) {
                return makeErrorMessage(ErrorCode::InvalidDemand, "Positive self-demand is not supported.");
            }

            total_demand += value;

            if (!std::isfinite(total_demand)) {
                return makeErrorMessage(ErrorCode::InvalidDemand, "Total demand overflow in link exposure analysis.");
            }
        }

        LinkExposureSummary result;

        const int edge_bound = graph.globalEdgeCount() * 2;

        for (int e = 0; e < edge_bound; e += 2) {
            if (!graph.edgeAlive(e))
                continue;

            const auto [u, v] = graph.getEdgeEndpoints(e);

            LinkExposureData link;
            link.edge_id = e;
            link.source = u;
            link.target = v;

            for (std::size_t i = 0; i < demand.size(); ++i) {
                const double value = demand.getDemandValue(i);

                if (value == 0.0)
                    continue;

                const auto [s, t] = demand.getDemandPair(i);

                const double fraction = scheme.getPhysicalLinkFlow(e, s, t);

                if (!std::isfinite(fraction) || fraction < 0.0) {
                    return makeErrorMessage(ErrorCode::InvalidRouting, "Invalid physical-link flow coefficient.");
                }

                if (fraction <= flow_epsilon)
                    continue;

                const double contribution = value * fraction;

                if (!std::isfinite(contribution)) {
                    return makeErrorMessage(ErrorCode::InvalidRouting, "Overflow while calculating link exposure.");
                }

                link.affected_demand_volume += value;
                link.lost_traffic += contribution;

                if (!std::isfinite(link.affected_demand_volume) || !std::isfinite(link.lost_traffic)) {
                    return makeErrorMessage(ErrorCode::InvalidRouting, "Link exposure accumulation overflow.");
                }
            }

            if (total_demand > 0.0) {
                link.lost_traffic_fraction = link.lost_traffic / total_demand;

                const double affected_fraction = link.affected_demand_volume / total_demand;

                result.maximum_affected_demand_fraction = std::max(result.maximum_affected_demand_fraction, affected_fraction);
            }

            if (link.lost_traffic > 0.0)
                ++result.traffic_carrying_links;

            if (link.lost_traffic_fraction >= 0.10)
                ++result.critical_links_10_percent;

            if (link.lost_traffic_fraction >= 0.25)
                ++result.critical_links_25_percent;

            if (link.lost_traffic_fraction >= 0.50)
                ++result.critical_links_50_percent;

            if (link.lost_traffic_fraction > result.maximum_lost_traffic_fraction || result.most_critical_edge_id == -1) {

                result.maximum_lost_traffic_fraction = link.lost_traffic_fraction;

                result.most_critical_edge_id = e;
                result.most_critical_source = u;
                result.most_critical_target = v;
            }

            result.links.push_back(link);
        }

        result.tested_links = result.links.size();

        if (result.links.empty())
            return result;

        std::vector<double> fractions;
        fractions.reserve(result.links.size());

        double sum = 0.0;

        for (const auto& link : result.links) {
            fractions.push_back(link.lost_traffic_fraction);
            sum += link.lost_traffic_fraction;

            if (!std::isfinite(sum)) {
                return makeErrorMessage(ErrorCode::InvalidRouting, "Exposure summary overflow.");
            }
        }

        result.average_lost_traffic_fraction = sum / static_cast<double>(fractions.size());

        std::sort(fractions.begin(), fractions.end());

        const std::size_t n = fractions.size();

        if (n % 2 == 1) {
            result.median_lost_traffic_fraction = fractions[n / 2];
        } else {
            // Avoid intermediate overflow from adding two doubles.
            result.median_lost_traffic_fraction = fractions[n / 2 - 1] / 2.0 + fractions[n / 2] / 2.0;
        }

        return result;
    }
};

#endif // E_ROUTING_LINK_EXPOSURE_ANALYZER_H
