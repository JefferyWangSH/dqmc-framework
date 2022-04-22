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
            void set_lattice_params(const SideLengthVec& side_length_vec);            

            // initializations
            void initial();

            // definition of vector products 
            // const double product(const std::array<double,2>& vecr, const std::array<double,2>& vecp);
        

        private:

            // private initialization functions
            void initial_hopping_matrix();
            void initial_index2site_table();
            void initial_nearest_neighbour_table();

    };
    
} // namespace Lattice

#endif // LATTICE_SQUARE_H