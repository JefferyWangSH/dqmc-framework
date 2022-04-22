#include "lattice/square.h"


namespace Lattice {

    double Square::product(const std::array<double,2>& vecr, const std::array<double,2>& vecp) {
        const auto& [rx, ry] = vecr;
        const auto& [px, py] = vecp;
        // for square lattice, the two basis are orthogonal
        return rx * px + ry * py;
    }

} // namespace Lattice