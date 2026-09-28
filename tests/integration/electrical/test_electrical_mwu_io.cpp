#include "../../common/utils.h"
#include <catch2/catch_approx.hpp>
#include <filesystem>


using namespace integration;

TEST_CASE("Electrical flow solver solves a tiny LGF dataset", "[integration][electrical][dataset][lgf]")
{
    const std::filesystem::path dataset = tinyLgfDataset();

    REQUIRE(std::filesystem::exists(dataset));

    Config cfg = makeElectricalConfig();
    cfg.filename = dataset;
    auto graph = makegraph(cfg.graph_format);
    if (!cfg.filename.empty()
        && graph) {
        auto r = GraphIO::readLGFFile(*(graph.value()), cfg.filename);
        if (!r) { FAIL(r.error().message.c_str()); }
    }

    (graph.value())->finalize();

    REQUIRE((graph.value()));
    REQUIRE((graph.value())->getNumNodes() > 0);
    REQUIRE((graph.value())->getNumUndirectedEdges() > 0);


    RoutingEngine engine;
    auto result = engine.solve(*((graph.value())), cfg, cfg.solvers.front());

    requireValidRoutingResult(result);
}

TEST_CASE("Electrical flow solver works with gravity demand model",
          "[integration][electrical][demand][gravity]")
{
    const std::filesystem::path dataset = tinyLgfDataset();

    REQUIRE(std::filesystem::exists(dataset));

    Config cfg = makeElectricalConfig();
    cfg.filename = dataset;
    auto graph = makegraph(cfg.graph_format);
    if (!cfg.filename.empty()
        && (graph)) {
        auto r = GraphIO::readLGFFile(*(graph.value()), cfg.filename);
        if (!r) { FAIL(r.error().message.c_str()); }
    }
    (graph.value())->finalize();

    cfg.demand_models.push_back(DemandModelType::GRAVITY); // adapt to your real enum/name
    cfg.evaluate_demand_models = true;
    cfg.offline_opt_per_model["gravity"]=computeOfflineOptimalCongestion(*(graph.value()), cfg.demand_maps["gravity"]);

    RoutingEngine engine;
    IRoutingResult result;
    if (auto res = engine.solve(*(graph.value()), cfg, cfg.solvers.front())) {
        result = std::move(res.value());
    }

    REQUIRE(std::isfinite(result.oblivious_ratio));
    REQUIRE(result.oblivious_ratio >= 1.0);
}

int runCommand(const std::string& command)
{
    return std::system(command.c_str());
}


TEST_CASE("Electrical solver CLI runs on a small LGF file with gravity demand",
          "[integration][cli][electrical]")
{
    const std::filesystem::path executable =
    std::filesystem::path(PROJECT_BINARY_DIR) / "oblivious_routing";

    const std::filesystem::path dataset = tinyLgfDataset();

    REQUIRE(std::filesystem::exists(executable));
    REQUIRE(std::filesystem::exists(dataset));

    const std::string command =
        executable.string() +
        " electrical " +
        dataset.string() +
        " gravity";

    const int exit_code = runCommand(command);

    REQUIRE(exit_code == 0);
}

TEST_CASE("Electrical solver CLI rejects missing demand model without crashing",
          "[integration][cli][electrical][error-handling]")
{
    const std::filesystem::path executable =
    std::filesystem::path(PROJECT_BINARY_DIR) / "oblivious_routing";

    const std::filesystem::path dataset = tinyLgfDataset();

    REQUIRE(std::filesystem::exists(executable));
    REQUIRE(std::filesystem::exists(dataset));

    const std::string command =
        executable.string() +
        " electrical " +
        dataset.string();

    const int exit_code = runCommand(command);

    REQUIRE(exit_code == 0);
}

TEST_CASE("Electrical solver conserves unit flow per source",
          "[integration][electrical][error-handling]") {
    const std::filesystem::path dataset = Backbone_1239_LgfDataset();

    REQUIRE(std::filesystem::exists(dataset));

    Config cfg = makeElectricalConfig();
    cfg.filename = dataset;
    auto graph = makegraph(cfg.graph_format);
    if (!cfg.filename.empty()
        && graph) {
        auto r = GraphIO::readLGFFile(*(graph.value()), cfg.filename);
        if (!r) { FAIL(r.error().message.c_str()); }
    }

    (graph.value())->finalize();

    ElectricalMWU solver(*(graph.value()), 0, true);
    std::unique_ptr<RoutingScheme> table;
    if (auto res = solver.solve()) {
        table = std::move(res.value());
    }

    constexpr double tol = EPS;

    for (int s = 0; s < (graph.value())->getNumNodes(); ++s) {
        double out = 0.0;
        double in = 0.0;
        for (int e = 0; e < (graph.value())->getNumDirectedEdges(); ++e) {
            auto [u, v] = (graph.value())->getEdgeEndpoints(e);
            double f = table->getFlow(e, s, 0);

            REQUIRE(std::isfinite(f));

            if (u == s) out += f;
            if (v == s) in += f;
        }

        REQUIRE((std::abs(out)-std::abs(in)) <= SOFT_EPS);
    }
}