#pragma once

#include "../hierarchy_results.h"

#include "core/errors.h"
#include "data_structures/graph/Igraph.h"

class XCutHierarchyPreprocessor {
public:

    Result<HierarchyResult> build(const IGraph& graph) const;

private:

    static std::pair<std::vector<std::pair<unsigned int,unsigned int>>, std::vector<double>>toXCutEdges(const IGraph& graph);

    static std::vector<int> computeInducedEdges(const IGraph& graph,const std::vector<int>& vertices);

    static void normalizeLevelOrder(HierarchyResult& hierarchy);

    static Result<void> buildParentChildRelations(HierarchyResult& hierarchy);

    static void buildLookupStructures(const IGraph& graph,HierarchyResult& hierarchy);

    static void choosePortals(
        const IGraph& graph,
        HierarchyResult& hierarchy
    );
};