//
// Created by Mert on 29.10.26.
//

#ifndef E_ROUTING_LINK_EXPOSURE_RESULT_ADAPTER_H
#define E_ROUTING_LINK_EXPOSURE_RESULT_ADAPTER_H
#include "analysis/failures/link_exposure_analyzer.h"
#include "routing/routing_result.h"

inline void populateLinkExposureMetrics(DemandEvaluationResult& evaluation, const LinkExposureSummary& exposure, double runtime_microseconds) {
    evaluation.failure_analysis_runtime_microseconds = runtime_microseconds;

    evaluation.failure_tested_links = exposure.tested_links;

    evaluation.failure_most_critical_edge_id = exposure.most_critical_edge_id;

    evaluation.failure_most_critical_source = exposure.most_critical_source;

    evaluation.failure_most_critical_target = exposure.most_critical_target;

    evaluation.failure_maximum_lost_traffic_fraction = exposure.maximum_lost_traffic_fraction;

    evaluation.failure_average_lost_traffic_fraction = exposure.average_lost_traffic_fraction;

    evaluation.failure_median_lost_traffic_fraction = exposure.median_lost_traffic_fraction;

    evaluation.failure_maximum_affected_demand_fraction = exposure.maximum_affected_demand_fraction;

    evaluation.failure_traffic_carrying_links = exposure.traffic_carrying_links;

    evaluation.failure_critical_links_10_percent = exposure.critical_links_10_percent;

    evaluation.failure_critical_links_25_percent = exposure.critical_links_25_percent;

    evaluation.failure_critical_links_50_percent = exposure.critical_links_50_percent;
}
#endif // E_ROUTING_LINK_EXPOSURE_RESULT_ADAPTER_H
