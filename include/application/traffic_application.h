//
// Created by Mert on 02.11.26.
//

#ifndef E_ROUTING_TRAFFIC_APPLICATION_H
#define E_ROUTING_TRAFFIC_APPLICATION_H
#include "io/cli/command.h"

class TrafficApplication {
    public:
    // Executes valid analyzed-traffic demand
    static int run(const Command& command);
};

#endif // E_ROUTING_TRAFFIC_APPLICATION_H
