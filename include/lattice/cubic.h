#ifndef LATTICE_CUBIC_H
#define LATTICE_CUBIC_H
#pragma once


/**
  *  This header file defines the Lattice::Cubic class for 3d cubic lattice,
  *  derived from the base class Lattice::LatticeBase.
  */

#include "lattice/lattice_base.h"

namespace Lattice {

    
    // ------------------------ Derived class Lattice::Cubic for 3d cubic lattice ----------------------------
    class Cubic : public LatticeBase {

        private:
            
            // // some high symmetric points in the reciprocal lattice
            // LatticeInt m_gamma_point_index{};
            // LatticeInt m_m_point_index{};
            // LatticeInt m_x_point_index{};
            // LatticeIntVec m_delta_line_index{};
            // LatticeIntVec m_z_line_index{};
            // LatticeIntVec m_sigma_line_index{};
            // LatticeIntVec m_gamma2x2m2gamma_loop_index{};   // indexes of the defined loop


        public:
        
            Cubic() = default;

            // set up lattice parameters
            void set_lattice_params(const LatticeIntVec& side_length_vec);            

            // initializations
            void initial();

            // // interfaces for high symmetric momentum points
            // const LatticeInt GammaPointIndex()     const ;
            // const LatticeInt MPointIndex()         const ;
            // const LatticeInt XPointIndex()         const ;
            // const LatticeIntVec& DeltaLineIndex()  const ;
            // const LatticeIntVec& ZLineIndex()      const ;
            // const LatticeIntVec& SigmaLineIndex()  const ;
            // const LatticeIntVec& Gamma2X2M2GammaLoopIndex() const ;
            

        private:

            // private initialization functions
            void initial_index2site_table();
            void initial_index2momentum_table();

            void initial_nearest_neighbour_table();
            void initial_displacement_table();
            void initial_symmetric_points();
            void initial_fourier_factor_table();

            void initial_hopping_matrix();

    };
    
} // namespace Lattice


#endif // LATTICE_CUBIC_H
