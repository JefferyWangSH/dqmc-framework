#ifndef MODEL_BASE_H
#define MODEL_BASE_H
#pragma once


/**
  *  This header file defines the abstract base class Model::ModelBase, 
  *  which is pure virtual, for describing different kinds of quantum models.
  *  Up to now, it only supports fermionic systems with spin-1/2,
  *  and the fermionic fields should interact with Z2 bosonic fields 
  *  after a real-valued Hubbard–Stratonovich transformation.
  */  

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


// forward declaration
namespace QuantumMonteCarlo {
    class DqmcWalker;
}
namespace Lattice {
    class LatticeBase;
}

namespace Model {

    // -------------------------------- Abstract base class Model::ModelBase ------------------------------------
    class ModelBase {
        protected:
            using SpaceIndex = int;
            using TimeIndex = int;
            using Spin = int;
            using Walker = QuantumMonteCarlo::DqmcWalker;
            using Lattice = Lattice::LatticeBase;
            using GreensFunc = Eigen::MatrixXd;
            using Matrix = Eigen::MatrixXd;
            
            // Exponential form of K matrix and V matrix
            // K matrix, corresponding to the hopping matrix, depending only on the hopping constants t's and the geometry of the lattice
            // V matrix corresponds to the interaction matrix between fermions and the auxiliary bosonic fields,
            // which emerges during the process of Hubbard–Stratonovich transformation in DQMC.
            // In principle V matrix should be model-dependent.

            // expK and expV(t) are defined as :
            //      expK = exp( -dt K )    and    expV(t) = exp(- dt V_sigma(t) )
            // where dt is the interval of the imaginary-time grids, 
            // and their transpose and inversion are straightforward.
            // cation that the V matrix depends on the configurations of bosonic fields and the spin of fermion,
            // and in general, it also depends on the model and specific HS transformation.
            Matrix m_expK_mat{};
            Matrix m_expV_mat{};
            Matrix m_inv_expK_mat{};
            Matrix m_inv_expV_mat{};
            Matrix m_trans_expK_mat{};
            Matrix m_trans_expV_mat{};

            // other model parameters should be defined in the derived Model classes .

        public:
            ModelBase() = default;
        
        protected:

            // initialize the model class for specific lattice
            virtual void initial(const Lattice& lattice, const Walker& walker) = 0;

            // return the updating radio of one step of the local dqmc update
            // which is model-dependent
            virtual double get_update_radio(const Walker& walker, SpaceIndex space_index, TimeIndex time_index) = 0;

            // perform one local dqmc update
            // it should be defined in this subroutine that 
            // how the equal-time green's functions are transformed 
            // given some specific updates of the bosonic fields
            virtual void local_update(const Walker& walker, SpaceIndex space_index, TimeIndex time_index) = 0;


            // B(t) matrix is defined as exp(- dt V_sigma(t) ) * exp( -dt K )
            // the following functions define the multiplications between green's functions and B matrices
            // which are frequently called when we wrap the green's functions to differnet imaginary-time slices.
            // from the perspective of effiency, we donot need to calculate expV explicitly,
            // because in many cases, depending on the specific model and HS transformation,
            // they are diagonalized and it is straightforward to directly define their muplication rules.
            virtual void mult_B_from_left       ( GreensFunc& green, const Walker& walker, TimeIndex time_index, Spin spin ) = 0;
            virtual void mult_B_from_right      ( GreensFunc& green, const Walker& walker, TimeIndex time_index, Spin spin ) = 0;
            virtual void mult_invB_from_left    ( GreensFunc& green, const Walker& walker, TimeIndex time_index, Spin spin ) = 0;
            virtual void mult_invB_from_right   ( GreensFunc& green, const Walker& walker, TimeIndex time_index, Spin spin ) = 0;
            virtual void mult_transB_from_left  ( GreensFunc& green, const Walker& walker, TimeIndex time_index, Spin spin ) = 0;
    };

}

#endif // MODEL_BASE_H