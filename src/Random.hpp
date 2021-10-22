#include <random>

/** global variable for random engine */
namespace Random {
    // random engine
    static std::default_random_engine Engine(time(nullptr));
}