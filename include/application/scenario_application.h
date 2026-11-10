//
// Created by Mert on 02.11.26.
//

#ifndef E_ROUTING_SCENARIO_DISPATCHER_H
#define E_ROUTING_SCENARIO_DISPATCHER_H

#include "io/cli/command.h"

class ScenarioApplication {
public:
    // Executes a validated analyze-scenario command.
    // Returns 0 on success and a nonzero code on failure.
    static int run(const Command& command);
};


#endif // E_ROUTING_SCENARIO_DISPATCHER_H
