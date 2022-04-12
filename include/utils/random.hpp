#ifndef UTIL_RANDOM_H
#define UTIL_RANDOM_H
#pragma once

#include <random>

namespace Utils {

    // -------------------- Utils::Random class for generating random seed used in MPI program ----------------------
    class Random {
        public:
            static std::default_random_engine engine;
            static void set_seed(const int &rank) {
                engine.seed( time(nullptr)+rank );
            }
    };

    // initialization of static member
    std::default_random_engine Random::engine( time(nullptr) );

} // namespace Utils

#endif // UTIL_RANDOM_H