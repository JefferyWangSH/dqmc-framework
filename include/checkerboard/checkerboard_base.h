#ifndef CHECKERBOARD_BASE_H
#define CHECKERBOARD_BASE_H
#pragma once


/**
  *  This header file defines the base class CheckerBoard::Base,
  *  for implementing the checkerboard breakups on biparticle lattice.
  *  It is pure virtual, and the breakups for any specific lattice should be derived from this class.
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace Model { class ModelBase; }
namespace Lattice { class LatticeBase; }
namespace QuantumMonteCarlo { class DqmcWalker; }

namespace CheckerBoard {

    using Site = std::vector<int>;
    using Matrix = Eigen::MatrixXd;
    using RealScalar = double;

    using ModelBase = Model::ModelBase;
    using LatticeBase = Lattice::LatticeBase;
    using DqmcWalker = QuantumMonteCarlo::DqmcWalker;


    // ------------------------------ Pure virtual base class CheckerBoard::Base --------------------------------
    class Base {
        public:
            
            // initialize from lattice, model and dqmcWalker
            virtual void set_params( const LatticeBase& lattice, 
                                     const ModelBase& model, 
                                     const DqmcWalker& walker ) = 0;
            virtual void initial() = 0;

            // multiply the exponent of hopping matrix K using checkerboard breakups
            virtual void mult_expK_from_left        ( Matrix &matrix ) const = 0;
            virtual void mult_expK_from_right       ( Matrix &matrix ) const = 0;
            virtual void mult_inv_expK_from_left    ( Matrix &matrix ) const = 0;
            virtual void mult_inv_expK_from_right   ( Matrix &matrix ) const = 0;
            virtual void mult_trans_expK_from_left  ( Matrix &matrix ) const = 0;

        protected:

            // multiply hopping matrix K within single plaquette, labeled by site vector
            virtual void mult_expK_plaquette_from_left       ( Matrix &matrix, const Site& site ) const = 0;
            virtual void mult_expK_plaquette_from_right      ( Matrix &matrix, const Site& site ) const = 0;
            virtual void mult_inv_expK_plaquette_from_left   ( Matrix &matrix, const Site& site ) const = 0;
            virtual void mult_inv_expK_plaquette_from_right  ( Matrix &matrix, const Site& site ) const = 0;

    };

} // namespace CheckerBoard 



#endif // CHECKERBOARD_BASE_H