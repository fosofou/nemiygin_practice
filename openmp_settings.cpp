#include "openmp_settings.h"
#include "omp.h"

OpenMPSettings::OpenMPSettings() : numThreads(omp_get_max_threads()) {}

int OpenMPSettings::getNumThreads() const {
    return numThreads;
}

void OpenMPSettings::setNumThreads(int threads) {
    if (threads > 0) {
        numThreads = threads;
        // omp_set_num_threads(threads);
    }
}
