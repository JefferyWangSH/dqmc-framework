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

        private:
            
            // some high symmetry points in the reciprocal lattice
            LatticeInt m_gamma_point_index{};           // (0,0)
            LatticeInt m_x_point_index{};               // (pi,0)
            LatticeInt m_m_point_index{};               // (pi,pi)
            LatticeIntVec m_delta_line_index{};         // (0,0)  ->  (pi,0)
            LatticeIntVec m_z_line_index{};             // (pi,0) ->  (pi,pi)
            LatticeIntVec m_sigma_line_index{};         // (0,0)  ->  (pi,pi)

            // defined loop (0,0) -> (pi,0) -> (pi,pi) -> (0,0) 
            LatticeIntVec m_gamma2x2m2gamma_loop_index{};   


        public:
        
            Square() = default;

            // set up lattice parameters
            void set_lattice_params(const LatticeIntVec& side_length_vec);            

            // initializations
            void initial();

            // interfaces for high symmetry momentum points
            const LatticeInt GammaPointIndex()     const ;
            const LatticeInt XPointIndex()         const ;
            const LatticeInt MPointIndex()         const ;
            const LatticeIntVec& DeltaLineIndex()  const ;
            const LatticeIntVec& ZLineIndex()      const ;
            const LatticeIntVec& SigmaLineIndex()  const ;
            const LatticeIntVec& Gamma2X2M2GammaLoopIndex() const ;
            

        private:

            // private initialization functions
            void initial_index2site_table();
            void initial_index2momentum_table();

            void initial_nearest_neighbour_table();
            void initial_displacement_table();
            void initial_symmetry_points();
            void initial_fourier_factor_table();

            void initial_hopping_matrix();

    };
    
} // namespace Lattice

#endif // LATTICE_SQUARE_H