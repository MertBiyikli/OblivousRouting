//
// Created by Mert Biyikli on 10.06.26.
//

#ifndef OBLIVIOUSROUTING_OPENMP_CHUNKED_H
#define OBLIVIOUSROUTING_OPENMP_CHUNKED_H

#include "../mwu/par_electrical_flow.h"

struct OpenMPChunkedExecution {
    int threads=1;

    int chunkSize = 4;

    int numWorkers() const {
        return threads;
    }

    template <typename Task>
    void parallelFor(int begin, int end, Task&& task) const {
#ifdef OR_ENABLE_OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic, chunkSize)
        for (int i = begin; i < end; ++i) {
            task(i, omp_get_thread_num());
        }
#else
        for (int i = begin; i < end; ++i) {
            task(i, 0);
        }
#endif
    }
};

#endif //OBLIVIOUSROUTING_OPENMP_CHUNKED_H