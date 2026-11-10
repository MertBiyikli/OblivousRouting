//
// Created by Mert on 30.10.26.
//

#include "io/cli/help.h"
#include "io/cli/options.h"
#include <sstream>

std::string solverList() {
    std::ostringstream out;

    out << "Available solvers:\n"
        << "\n"
        << "  electrical       "
        << "Electrical Flow (sketching)\n"

        << "  electrical_naive "
        << "Electrical Flow (naive/exact loads)\n"

        << "  raecke_frt       "
        << "Raecke MWU + FRT\n"

        << "  raecke_ckr       "
        << "Raecke MWU + Fast-CKR\n"

        << "  raecke_mst       "
        << "Raecke MWU + randomized MST\n"

        << "  cohen            "
        << "Applegate-Cohen LP\n"

        << "  expander         "
        << "Electrified Expander Hierarchy\n"

        << "  expander_mwu     "
        << "Expander hierarchy + MWU\n"

        << "  semi_electrical  "
        << "Semi-oblivious routing, electrical base\n"

        << "  semi_tree        "
        << "Semi-oblivious routing, tree base\n"

        << "  semi_expander    "
        << "Semi-oblivious routing, expander base\n";

    return out.str();
}

std::string help(std::string_view executable) {
    std::ostringstream out;

    const auto options = solveOptions();

    out << "E-Routing: Network Resilience and Capacity Planning Engine.\n"
        << "\n"

        << "Usage:\n"

        << "  " << executable << " solve "
        << "--solver <name[,name...]> "
        << "--graph <file> "
        << "[options]\n"

        << "  " << executable << " analyze-n1 --solver <name> --graph <file> --demand <model> --output <report.json> [options]\n"

        << "  " << executable << " --list-solvers\n"

        << "  " << executable << " --help\n"

        << "\n"

        << options

        << "\n"

        << "Examples:\n"

        << "  " << executable << " solve "
        << "--solver electrical "
        << "--graph "
        << "experiments/datasets/small/Backbone/1221.lgf"
        << "\n"

        << "\n"

        << "  " << executable << " solve "
        << "-s electrical,raecke_ckr "
        << "-g experiments/datasets/small/Backbone/1221.lgf "
        << "-d uniform "
        << "-j 4 "
        << "-o result/backbone.json"
        << "\n";

    return out.str();
}
