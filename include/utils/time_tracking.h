//
// Created by Mert Biyikli on 20.03.26.
//

#ifndef OBLIVIOUSROUTING_TIME_TRACKING_H
#define OBLIVIOUSROUTING_TIME_TRACKING_H

#include <chrono>

#define duration(a) std::chrono::duration_cast<std::chrono::microseconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

#endif //OBLIVIOUSROUTING_TIME_TRACKING_H
