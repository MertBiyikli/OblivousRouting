//
// Created by Mert Biyikli on 24.07.26.
//

#ifndef OBLIVIOUSROUTING_JSON_EXPORTER_H
#define OBLIVIOUSROUTING_JSON_EXPORTER_H

#include <filesystem>
#include <string>

#include "visualization_result.h"
#include "core/errors.h"

class RoutingVisualizationJsonExporter {
public:
    static Result<void> write(const RoutingVisualizationResult& result,const std::filesystem::path& output_path);
    static std::string serialize(const RoutingVisualizationResult& result);
};

#endif //OBLIVIOUSROUTING_JSON_EXPORTER_H