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

#include <memory>
#include <functional>
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

namespace CheckerBoard {
    class CheckerBoardBase;
}

namespace Model {

    // useful aliases
    using SpaceIndex = int;
    using TimeIndex = int;
    using Spin = int;
    using RealScalar = double;

    using Walker = QuantumMonteCarlo::DqmcWalker;
    using LatticeBase = Lattice::LatticeBase;
    using CheckerBoardBase = CheckerBoard::CheckerBoardBase;
    using Matrix = Eigen::MatrixXd;

    using GreensFunc = Eigen::MatrixXd;
    using GreensFuncVec = std::vector<Eigen::MatrixXd>;
    using ptrGreensFunc = std::unique_ptr<Eigen::MatrixXd>;
    using ptrGreensFuncVec = std::unique_ptr<std::vector<Eigen::MatrixXd>>;

    using MultExpKMethod = void( GreensFunc& );
    using MultInvExpKMethod = void( GreensFunc& );
    using MultTransExpKMethod = void( GreensFunc& );


    // ------------------------------------- Abstract base class Model::ModelBase ------------------------------------------
    class ModelBase {
        protected:

            // size of the space-time lattices
            int m_space_size{};
            int m_time_size{};

            // K matrix, corresponding to the hopping matrix, 
            // depends only on the hopping constants t's and the geometry of the lattice.
            // V matrix corresponds to the interaction matrix between fermions and the auxiliary bosonic fields,
            // which emerges during the process of Hubbard–Stratonovich transformation in DQMC.
            // In principle V matrix should be model-dependent.

            // expK and expV(t) are defined as :
            //      expK = exp( -dt K )    and    expV(t) = exp(- dt V_sigma(t) )
            // where dt is the interval of the imaginary-time grids, 
            // and their transpose and inversion are straightforward.
            // cation that the V matrix depends on the configurations of bosonic fields and the spin of fermion,
            // and in general, it also depends on the model and specific HS transformation.

            // From the perspective of effiency, we donot need to calculate expV explicitly,
            // because in many cases, depending on the specific model and HS transformation,
            // they are diagonalized and it is straightforward to directly define their muplication rules.
            // Meanwhile, for specific quantum model on a biparticle lattice, 
            // it should be much more efficient if checkerboard breakups can be applied,
            // and in that case there is also no need to explicitly compute the exponential of K.
            
            // In conclusion, directly computing the expK or expV are not a necessity
            // For better performance, depending on specific lattice and model, 
            // it's recommended to link the Model class to a specific checkerboard class.
            Matrix m_expK_mat{};
            Matrix m_expV_mat{};
            Matrix m_inv_expK_mat{};
            Matrix m_inv_expV_mat{};
            Matrix m_trans_expK_mat{};
            Matrix m_trans_expV_mat{};

            // function pointers for multiplying exponent of K matrix to a dense matrix
            // these cound be redirected to any checkerboard methods 
            // to get accelerated by checkerboard breakups
            std::function<MultExpKMethod>       m_mult_expK_from_left{};
            std::function<MultExpKMethod>       m_mult_expK_from_right{};
            std::function<MultInvExpKMethod>    m_mult_inv_expK_from_left{};
            std::function<MultInvExpKMethod>    m_mult_inv_expK_from_right{};
            std::function<MultTransExpKMethod>  m_mult_trans_expK_from_left{};

            // other model parameters should be defined in the derived Model classes .


        public:

            ModelBase() = default;

            // ------------------------------------------ Setup interfaces -----------------------------------------------

            // setup model params ( interface for derived classed )
            virtual void set_model_params(RealScalar, RealScalar, RealScalar) = 0;

            virtual const RealScalar HoppingT() const = 0; 
            virtual const RealScalar ChemicalPotential() const = 0;


            // ----------------------------------------- Initializations -------------------------------------------------

            // initialize the model class for specific lattice and DqmcWalker
            virtual void initial             (const LatticeBase& lattice, const Walker& walker) = 0;
            virtual void initial_params      (const LatticeBase& lattice, const Walker& walker) = 0;
            virtual void initial_KV_matrices (const LatticeBase& lattice, const Walker& walker) = 0;

            // randomrize the bosonic fields, which is model-dependent
            virtual void set_bosonic_fields_to_random() = 0;

            
            // ------------------------------------------ Linking methods ------------------------------------------------

            // naively link member functions mult_expK's to m_mult_expK methods
            // this will multply the exponent of K directly to a dense matrix
            void link();

            // link checkerboard member functions to m_mult_expK methods
            // this will achieve higher performance by using checkerboard breakups on a specific lattice
            void link( const CheckerBoardBase& checkerboard );


            // ---------------------------------------- Monte Carlo updates ----------------------------------------------

            // perform one local dqmc update
            virtual void update_bosonic_field(TimeIndex time_index, SpaceIndex space_index) = 0;

            // return the updating ratio of one step of the local dqmc update
            // which is model-dependent
            virtual const double get_update_ratio(Walker& walker, TimeIndex time_index, SpaceIndex space_index) const = 0;
            
            // transform the equal-time green's functions
            // given a specific update of the bosonic fields
            virtual void update_greens_function(Walker& walker, TimeIndex time_index, SpaceIndex space_index) = 0;
            

            // ------------------------------------------ Warpping methods -----------------------------------------------

            // B(t) matrix is defined as exp(- dt V_sigma(t) ) * exp( -dt K )
            // the following functions define the multiplications between green's functions and B matrices
            // which are frequently called when we wrap the green's functions to differnet imaginary-time slices.
            // Caution: these are perfermance-important steps.
            virtual void mult_B_from_left       ( GreensFunc& green, TimeIndex time_index, Spin spin ) const = 0;
            virtual void mult_B_from_right      ( GreensFunc& green, TimeIndex time_index, Spin spin ) const = 0;
            virtual void mult_invB_from_left    ( GreensFunc& green, TimeIndex time_index, Spin spin ) const = 0;
            virtual void mult_invB_from_right   ( GreensFunc& green, TimeIndex time_index, Spin spin ) const = 0;
            virtual void mult_transB_from_left  ( GreensFunc& green, TimeIndex time_index, Spin spin ) const = 0;

        
        private:

            // naive implementation of multplying exponent of K matrices to the greens function
            // note that these member function should not be explicitly called
            // using std::function member instead.
            void mult_expK_from_left        ( GreensFunc& green ) const ;
            void mult_expK_from_right       ( GreensFunc& green ) const ;
            void mult_inv_expK_from_left    ( GreensFunc& green ) const ;
            void mult_inv_expK_from_right   ( GreensFunc& green ) const ;
            void mult_trans_expK_from_left  ( GreensFunc& green ) const ;
    };

}

#endif // MODEL_BASE_H