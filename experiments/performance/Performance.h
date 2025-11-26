//
// Created by Mert Biyikli on 12.08.25.
//

#ifndef OBLIVIOUSROUTING_PERFORMANCE_H
#define OBLIVIOUSROUTING_PERFORMANCE_H

#include <chrono>
#include <functional>
#include "../../src/datastructures/graph.h"

class Performance{
public:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = Clock::time_point;

    static TimePoint start() {
        return Clock::now();
    }

    static double elapsed(TimePoint start) {
        auto end = Clock::now();
        return std::chrono::duration<double>(end - start).count();
    }

    static void printElapsed(const std::string& message, TimePoint start) {
        double duration = elapsed(start);
        std::cout << message << " took " << duration << " seconds." << std::endl;
    }

    void run(std::function<void()> func, Graph& g, std::string& message) {
        TimePoint start = this->start();
        func();
        printElapsed(message, start);
    }
};

#endif //OBLIVIOUSROUTING_PERFORMANCE_H
