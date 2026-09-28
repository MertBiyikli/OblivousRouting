//
// Created by Mert Biyikli on 10.06.26.
//

#ifndef OBLIVIOUSROUTING_OPENMP_DYNAMIC_H
#define OBLIVIOUSROUTING_OPENMP_DYNAMIC_H

#include "../mwu/par_electrical_flow.h"

struct OpenMPDynamicExecution {
    int threads=1;


    int numWorkers() const {
        return threads;
    }

    template <typename Task>
    void parallelFor(int begin, int end, Task&& task) const {
#ifdef OR_ENABLE_OPENMP
        #pragma omp parallel num_threads(threads)
        {
            omp_set_dynamic(0);
            omp_set_max_active_levels(1);

            #pragma omp for schedule(dynamic, 1)
            for (int i = begin; i < end; ++i) {
                task(i, omp_get_thread_num());
            }
        }
#else
        for (int i = begin; i < end; ++i) {
            task(i, 0);
        }
#endif
    }
};

#endif //OBLIVIOUSROUTING_OPENMP_DYNAMIC_H