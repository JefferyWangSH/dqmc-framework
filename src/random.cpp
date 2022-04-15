#include "random.h"

namespace Utils {

    // initialization of static member
    std::default_random_engine Random::Engine( time(nullptr) );

    void Random::set_seed(const int& rank) {
        Engine.seed( time(nullptr)+rank );
    }

    void Random::set_seed_fixed(const int& seed) {
        Engine.seed( seed );
    }

} // namespace Utils
