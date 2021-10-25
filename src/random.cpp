#include "random.h"

namespace Random {
    // definition
    std::default_random_engine Engine(time(nullptr));

    // set up random seed for different processors
    void set_seed(const int &rank) {
        Engine.seed( time(nullptr)+rank );
    } 
}
