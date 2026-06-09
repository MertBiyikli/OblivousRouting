#include "../../mwu/electrical_mwu.h"

class ParElectricalFlowMWU : public ElectricalMWU {
    ParElectricalFlowMWU(IGraph& g, int root, bool use_sketching, bool debug = false)
    : ElectricalMWU(g, root, use_sketching, debug){
    }
};