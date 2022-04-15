#ifndef UTIL_RANDOM_H
#define UTIL_RANDOM_H
#pragma once

#include <random>

namespace Utils {

    // -------------------- Utils::Random class for generating random seed used in MPI program ----------------------
    class Random {
        public:
            static std::default_random_engine Engine;

            // setup seeds randomly accroding to the processor rank in MPI
            static void set_seed(const int& rank);

            // setup fixed seeds for debug usages
            static void set_seed_fixed(const int& seed);
    };

} // namespace Utils

#endif // UTIL_RANDOM_H