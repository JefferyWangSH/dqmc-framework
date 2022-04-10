#ifndef LATTICE_SQUARE_3D_H
#define LATTICE_SQUARE_3D_H
#pragma once

/**
  *  This header file defines the Lattice::Square3d class for 3d square lattice,
  *  inherited from the base class Lattice::LatticeBase.
  */


#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "lattice/lattice_base.h"

namespace Lattice {

    class Square3d : public LatticeBase {
        public:
            Square3d() = default;
            Square3d(int space_size);

            void initial();
    };
    
} // namespace Lattice

#endif // LATTICE_SQUARE_3D_H