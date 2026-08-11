#pragma once

#include "../hierarchy_results.h"

#include "core/errors.h"
#include "data_structures/graph/Igraph.h"
#include "data_structures/graph/graph.h"

class XCutHierarchyPreprocessor {
public:

    Result<HierarchyResult> build(const optimized::Graph<EdgeData>& graph) const;

private:

    static std::pair<std::vector<std::pair<unsigned int,unsigned int>>, std::vector<double>>toXCutEdges(const optimized::Graph<EdgeData>& graph);

    static std::vector<int> computeInducedEdges(const optimized::Graph<EdgeData>& graph,const std::vector<int>& vertices);

    static void normalizeLevelOrder(HierarchyResult& hierarchy);

    static Result<void> buildParentChildRelations(HierarchyResult& hierarchy);

    static void buildLookupStructures(const optimized::Graph<EdgeData>& graph,HierarchyResult& hierarchy);

    static void choosePortals(const optimized::Graph<EdgeData>& graph,HierarchyResult& hierarchy);
};