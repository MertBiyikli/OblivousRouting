//
// Created by Mert on 02.11.26.
//

#ifndef E_ROUTING_RESILIENCE_APPLICATION_H
#define E_ROUTING_RESILIENCE_APPLICATION_H


#include "io/cli/command.h"

class ResilienceApplication {
public:
    // Executes a validated analyze-n1 command.
    // Returns 0 on success and nonzero on failure.
    static int run(const Command& command);
};

#endif // E_ROUTING_RESILIENCE_APPLICATION_H
