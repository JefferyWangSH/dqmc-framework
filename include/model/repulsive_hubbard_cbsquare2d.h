#ifndef REPULSIVE_HUBBARD_CBSQUARE_2D_H
#define REPULSIVE_HUBBARD_CBSQUARE_2D_H
#pragma once

/**
  *  This header file defines the specialized model class Model::RepulsiveHubbardCbSquare2d,
  *  designed for the repulsive hubbard model on 2d square lattice using checkboard breakups.
  *  The B-matrix multiplication methods are reloaded to achieve higher performance.
  */
 

#include "model/repulsive_hubbard.h"
#include "checkerboard/square2d.h"

namespace CheckerBoard {
    class Square2d;
}

namespace Model {

    // ------------------------- Specialized derived class Model::RepulsiveHubbardCbSquare2d ---------------------------
    class RepulsiveHubbardCbSquare2d : public RepulsiveHubbard {

        public:

            // checkerboard class for 2d square lattice
            CheckerBoard::Square2d m_checkerboard{};

            // initialize checkerboard class
            void initial_checkerboard( const LatticeBase& lattice, const Walker& walker );


        public:

            void initial( const LatticeBase& lattice, const Walker& walker );

            void mult_B_from_left       ( GreensFunc& green, TimeIndex time_index, Spin spin ) const ;
            void mult_B_from_right      ( GreensFunc& green, TimeIndex time_index, Spin spin ) const ;
            void mult_invB_from_left    ( GreensFunc& green, TimeIndex time_index, Spin spin ) const ;
            void mult_invB_from_right   ( GreensFunc& green, TimeIndex time_index, Spin spin ) const ;
            void mult_transB_from_left  ( GreensFunc& green, TimeIndex time_index, Spin spin ) const ;
    };

}

#endif // REPULSIVE_HUBBARD_CBSQUARE_2D_H
