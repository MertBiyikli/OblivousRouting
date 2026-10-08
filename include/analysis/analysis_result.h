//
// Created by Mert on 08.10.26.
//

#ifndef E_ROUTING_ANALYSIS_RESULT_H
#define E_ROUTING_ANALYSIS_RESULT_H

#include <map>
#include <optional>
#include <string>
#include <vector>


/*
 * This file defines the data structures used to represent the analysis result of a network routing analysis.
 * The structures are designed to be serializable and can be used to store and retrieve analysis results in a structured format.
 */


/**
 * @brief Represents an error that occurred during the analysis.
 */
struct AnalysisError {
    std::string errorCode;
    std::string errorMessage;
    bool redoable = false;
};

/**
 * @brief Represents the evidence information of the analysis.
 */
struct Evidence {
    std::string result_path;
    std::string statement_type;
};

/**
 * @brief Represents the metadata information of the analysis.
 */
struct AnalysisMetadata {
    std::string run_id, scenario_id, created_at, status;
    std::vector<std::string> enabled_analyses;
    std::optional<AnalysisError> error;
};

/**
 * @brief Represents the topology information of the network.
 */
struct Topology {
    std::string topology_id, sha256, source;
    std::uint64_t node_count = 0, physical_link_count = 0;
};

/**
 * @brief Represents the demand information of the network.
 */
struct Demand {
    std::string demand_id, sha256, source_type;
    std::optional<std::string> model;
    double scale_factor = 1.0, total_demand = 0;
    std::string unit;
};

/**
 * @brief Represents the link metric information of the network.
 */
struct LinkMetric {
    std::string link_id;
    double load = 0, capacity = 0, utilization = 0, headroom = 0;
    bool overloaded = false;
};

/**
 * @brief Represents the baseline metrics of the network.
 */
struct BaselineMetrics {
    double max_congestion = 0, average_utilization = 0;
    std::vector<LinkMetric> links;
    std::vector<std::string> overloaded_link_ids;
    std::uint64_t solver_runtime_microseconds = 0, evaluation_runtime_microseconds = 0;
};

/**
 * @brief Represents the baseline information of the network.
 */
struct Baseline {
    std::string status;
    std::optional<BaselineMetrics> metrics;
};

/**
 * @brief Represents the capacity headroom information of the network.
 */
struct CapacityHeadroom {
    std::string status;
    std::vector<double> thresholds;
    std::optional<double> minimum_headroom;
    std::vector<std::string> violating_link_ids;
};
/**
 * @brief Represents the exposure information of the network links.
 */
struct ExposureLink {
    std::string link_id;
    double traffic_exposed = 0, traffic_fraction = 0, affected_demand_pair_fraction = 0;
};
/**
 * @brief Represents the layer 1 exposure information of the network.
 */
struct Layer1Exposure {
    std::string status;
    std::vector<ExposureLink> links;
};
/**
 * @brief Represents a failure event in the network.
 */
struct FailureEvent {
    std::string event_id, failed_link_id, status, connectivity;
    double total_demand = 0, unroutable_demand = 0;
    std::optional<double> baseline_congestion, post_failure_congestion, degradation_factor;
    std::vector<std::string> bottleneck_link_ids;
    bool recomputation_success = false;
    std::optional<std::uint64_t> recomputation_runtime_microseconds;
    std::vector<std::string> capacity_violation_link_ids, affected_flow_ids;
    std::optional<AnalysisError> error;
};
/**
 * @brief Represents the N-1 recovery information of the network.
 */
struct N1Recovery {
    std::string status;
    std::vector<FailureEvent> events;
    std::uint64_t tested_link_count = 0;
    std::optional<double> service_survivability_rate, capacity_compliance_rate;
};
/**
 * @brief Represents a critical link in the network.
 */
struct CriticalLink {
    std::string link_id, metric;
    double value = 0;
    std::uint64_t rank = 0;
    std::vector<Evidence> evidence;
};
/**
 * @brief Represents a capacity violation in the network.
 */
struct CapacityViolation {
    std::string scenario_id, link_id;
    double utilization = 0, threshold = 0;
};
/**
 * @brief Represents a named scenario in the network.
 */
struct NamedScenario {
    std::string scenario_id, status;
};
/**
 * @brief Represents a growth scenario in the network.
 */
struct GrowthScenario {
    std::string scenario_id;
    double demand_scale_factor = 1;
    std::string status;
};
/**
 * @brief Represents an upgrade candidate in the network.
 */
struct UpgradeCandidate {
    std::string scenario_id, description, status;
};
/**
 * @brief Represents a recommendation in the network.
 */
struct Recommendation {
    std::string text, classification;
    std::vector<Evidence> evidence;
};
/**
 * @brief Represents a tool call in the network.
 */
struct ToolCall {
    std::string operation, status;
    std::map<std::string, std::string> parameters;
    std::optional<std::string> result_path;
};
/**
 * @brief Represents the provenance information of the network.
 */
struct Provenance {
    std::string engine_version, engine_commit, solver_id, solver_version, determinism;
    std::uint64_t random_seed = 0, thread_count = 1;
    std::map<std::string, std::string> configuration;
    std::vector<ToolCall> tool_calls;
};
/**
 * @brief Represents the analysis result of the network.
 */
struct AnalysisResult {
    std::string schema_version = "1.0.0";
    AnalysisMetadata analysis_metadata;
    Topology topology;
    Demand demand;
    Baseline baseline;
    CapacityHeadroom capacity_headroom;
    Layer1Exposure layer1_exposure;
    N1Recovery n1_recovery;
    std::vector<CriticalLink> critical_links;
    std::vector<CapacityViolation> capacity_violations;
    std::vector<std::string> disconnecting_failures;
    std::vector<NamedScenario> maintenance_scenarios;
    std::vector<GrowthScenario> growth_scenarios;
    std::vector<UpgradeCandidate> upgrade_candidates;
    std::vector<Recommendation> recommendations;
    Provenance provenance;
};
#endif // E_ROUTING_ANALYSIS_RESULT_H
