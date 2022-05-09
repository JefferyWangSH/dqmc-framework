#ifndef LATTICE_SQUARE_H
#define LATTICE_SQUARE_H
#pragma once

/**
  *  This header file defines the Lattice::Square class for 2d square lattice,
  *  derived from the base class Lattice::LatticeBase.
  */

#include "lattice/lattice_base.h"

namespace Lattice {

    
    // ------------------------ Derived class Lattice::Square for 2d square lattice ----------------------------
    class Square : public LatticeBase {
        public:
        
            Square() = default;

            // set up lattice parameters
            void set_lattice_params(const LatticeIntVec& side_length_vec);            

            // initializations
            void initial();
        

        private:

            // private initialization functions
            void initial_hopping_matrix();
            void initial_index2site_table();
            void initial_nearest_neighbour_table();
            void initial_index2momentum_table();
            void initial_symmetric_points();
            void initial_fourier_factor_table();

    };
    
} // namespace Lattice

#endif // LATTICE_SQUARE_H