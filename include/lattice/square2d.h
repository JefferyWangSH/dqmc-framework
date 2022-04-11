#ifndef LATTICE_SQUARE_2D_H
#define LATTICE_SQUARE_2D_H
#pragma once

/**
  *  This header file defines the Lattice::Square2d class for 2d square lattice,
  *  inherited from the base class Lattice::LatticeBase.
  */


#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "lattice/lattice_base.h"

namespace Lattice {
    
    // --------------------Derived class Lattice::Square2d for 2d square lattice ----------------------------
    class Square2d : public LatticeBase {
        public:
            Square2d() = default;

            Square2d(int space_size) : LatticeBase(space_size){};

            double product(const std::array<double,2>& vecr, const std::array<double,2>& vecp);
    };
    
} // namespace Lattice

#endif // LATTICE_SQUARE_2D_H