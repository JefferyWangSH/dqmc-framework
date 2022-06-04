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
            
            // some high symmetry points in the reciprocal lattice
            LatticeInt m_gamma_point_index{};           // (0,0,0)
            LatticeInt m_x_point_index{};               // (pi,0,0)
            LatticeInt m_m_point_index{};               // (pi,pi,0)
            LatticeInt m_r_point_index{};               // (pi,pi,pi)
            LatticeIntVec m_delta_line_index{};         // (0,0,0)   ->  (pi,0,0)
            LatticeIntVec m_z_line_index{};             // (pi,0,0)  ->  (pi,pi,0)
            LatticeIntVec m_sigma_line_index{};         // (0,0,0)   ->  (pi,pi,0)
            LatticeIntVec m_lambda_line_index{};        // (0,0,0)   ->  (pi,pi,pi)
            LatticeIntVec m_s_line_index{};             // (pi,0,0)  ->  (pi,pi,pi)
            LatticeIntVec m_t_line_index{};             // (pi,pi,0) ->  (pi,pi,pi)


        public:
        
            Cubic() = default;

            // set up lattice parameters
            void set_lattice_params(const LatticeIntVec& side_length_vec);            

            // initializations
            void initial();

            // interfaces for high symmetry momentum points
            const LatticeInt GammaPointIndex()     const ;
            const LatticeInt XPointIndex()         const ;
            const LatticeInt MPointIndex()         const ;
            const LatticeInt RPointIndex()         const ;
            const LatticeIntVec& DeltaLineIndex()  const ;
            const LatticeIntVec& ZLineIndex()      const ;
            const LatticeIntVec& SigmaLineIndex()  const ;
            const LatticeIntVec& LambdaLineIndex() const ;
            const LatticeIntVec& SLineIndex()      const ;
            const LatticeIntVec& TLineIndex()      const ;
            

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


#endif // LATTICE_CUBIC_H
