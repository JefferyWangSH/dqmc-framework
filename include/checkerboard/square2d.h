#ifndef CHECKERBOARD_SQUARE_2D_H
#define CHECKERBOARD_SQUARE_2D_H
#pragma once


/**
  *  This header file defines the CheckerBoard::Square2d class for the checkerboard breakups
  *  of 2d square lattice, which is derived from the base class CheckerBoard::Base .
  *  Notice that the break-ups can only be applied to the square lattice with even side length.
  */

#include "checkerboard/checkerboard_base.h"


namespace CheckerBoard {

    class Square2d : public Base {
        private:

            int m_side_length{};                // side length of the lattice
            int m_space_size{};                 // total number of sites
            RealScalar m_hopping_t{};
            RealScalar m_chemical_potential{};
            RealScalar m_time_interval{};

            // exponent of (inverse) hopping within plaquette
            Eigen::Matrix4d m_expK_plaquette{};
            Eigen::Matrix4d m_inv_expK_plaquette{};


        public:
            // set up parameters
            void set_params( const LatticeBase& lattice, 
                             const ModelBase& model, 
                             const DqmcWalker& walker );

            // initialization
            void initial();

            // multiply the exponent of hopping matrix K using checkerboard breakups
            void mult_expK_from_left        ( Matrix &matrix ) const ;
            void mult_expK_from_right       ( Matrix &matrix ) const ;
            void mult_inv_expK_from_left    ( Matrix &matrix ) const ;
            void mult_inv_expK_from_right   ( Matrix &matrix ) const ;
            void mult_trans_expK_from_left  ( Matrix &matrix ) const ;


        private:
            // multiply hopping matrix K within single plaquette, labeled by site vector
            void mult_expK_plaquette_from_left       ( Matrix &matrix, const Site& site ) const ;
            void mult_expK_plaquette_from_right      ( Matrix &matrix, const Site& site ) const ;
            void mult_inv_expK_plaquette_from_left   ( Matrix &matrix, const Site& site ) const ;
            void mult_inv_expK_plaquette_from_right  ( Matrix &matrix, const Site& site ) const ;
    };


} // namespace CheckerBoard


#endif // CHECKERBOARD_SQUARE_2D_H