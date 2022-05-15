#include "random.h"

namespace Utils {

    // initialization of the static member
    std::default_random_engine Random::Engine( time(nullptr) );

    void Random::set_seed(const int seed) {
        Engine.seed( seed );
    }

} // namespace Utils
