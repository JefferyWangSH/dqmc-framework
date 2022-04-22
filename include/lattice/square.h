#ifndef LATTICE_SQUARE_H
#define LATTICE_SQUARE_H
#pragma once

/**
  *  This header file defines the Lattice::Square class for 2d square lattice,
  *  derived from the base class Lattice::LatticeBase.
  */

#include "lattice/lattice_base.h"

namespace Lattice {
    
    // --------------------Derived class Lattice::Square for 2d square lattice ----------------------------
    class Square : public LatticeBase {
        public:
            Square() = default;

            Square(int space_size) : LatticeBase(space_size){};

            double product(const std::array<double,2>& vecr, const std::array<double,2>& vecp);
    };
    
} // namespace Lattice

#endif // LATTICE_SQUARE_H