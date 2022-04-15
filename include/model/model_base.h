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

namespace Measure {
    class MeasureHandler;
}

namespace Model {

    // useful aliases
    using SpaceIndex = int;
    using TimeIndex = int;
    using Spin = int;

    using Walker = QuantumMonteCarlo::DqmcWalker;
    using Lattice = Lattice::LatticeBase;
    using MeasureHandler = Measure::MeasureHandler;
    using Matrix = Eigen::MatrixXd;

    using GreensFunc = Eigen::MatrixXd;
    using GreensFuncVec = std::vector<Eigen::MatrixXd>;
    using ptrGreensFunc = std::unique_ptr<Eigen::MatrixXd>;
    using ptrGreensFuncVec = std::unique_ptr<std::vector<Eigen::MatrixXd>>;


    // -------------------------------- Abstract base class Model::ModelBase ------------------------------------
    class ModelBase {
        protected:

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
            Matrix m_expK_mat{};
            Matrix m_expV_mat{};
            Matrix m_inv_expK_mat{};
            Matrix m_inv_expV_mat{};
            Matrix m_trans_expK_mat{};
            Matrix m_trans_expV_mat{};

            // Equal-time green's functions, which is the most crucial quantities during dqmc simulations
            // for spin-1/2 systems, we label the spin index with up and down
            ptrGreensFunc m_green_tt_up{}, m_green_tt_dn{};
            ptrGreensFuncVec m_vec_green_tt_up{}, m_vec_green_tt_dn{};

            // Time-displaced green's functions G(t,0) and G(0,t)
            // important for time-displaced measurements of physical observables
            ptrGreensFunc m_green_t0_up{}, m_green_t0_dn{};
            ptrGreensFunc m_green_0t_up{}, m_green_0t_dn{};
            ptrGreensFuncVec m_vec_green_t0_up{}, m_vec_green_t0_dn{};
            ptrGreensFuncVec m_vec_green_0t_up{}, m_vec_green_0t_dn{};

            // bool params to indicate whether the equal-time or dynamic measurements will be performed
            bool m_is_equaltime{};
            bool m_is_dynamic{};

            int m_space_size{};
            int m_time_size{};

            // other model parameters should be defined in the derived Model classes .

        public:
            ModelBase() = default;

            // interface for protected members
            GreensFunc& GreenttUp();
            GreensFunc& GreenttDn();
            GreensFunc& Greent0Up();
            GreensFunc& Greent0Dn();
            GreensFunc& Green0tUp();
            GreensFunc& Green0tDn();

            GreensFuncVec& vecGreenttUp();
            GreensFuncVec& vecGreenttDn();
            GreensFuncVec& vecGreent0Up();
            GreensFuncVec& vecGreent0Dn();
            GreensFuncVec& vecGreen0tUp();
            GreensFuncVec& vecGreen0tDn();
        
        protected:

            // initialize the model class for specific lattice and DqmcWalker
            virtual void initial_KV_matrices(const Lattice& lattice, const Walker& walker) = 0;
            void initial_greens_function(const Lattice& lattice, const Walker& walker, const MeasureHandler& meas_handler);
            virtual void initial(const Lattice& lattice, const Walker& walker, const MeasureHandler& meas_handler) = 0;

            // randomrize the bosonic fields, which is model-dependent
            virtual void set_bosonic_fields_to_random() = 0;

            // return the updating radio of one step of the local dqmc update
            // which is model-dependent
            virtual const double get_update_radio(TimeIndex time_index, SpaceIndex space_index) const = 0;

            // perform one local dqmc update
            virtual void update_bosonic_field(TimeIndex time_index, SpaceIndex space_index) = 0;
            
            // transform the equal-time green's functions
            // given a specific update of the bosonic fields
            virtual void update_greens_function(TimeIndex time_index, SpaceIndex space_index) = 0;

            // B(t) matrix is defined as exp(- dt V_sigma(t) ) * exp( -dt K )
            // the following functions define the multiplications between green's functions and B matrices
            // which are frequently called when we wrap the green's functions to differnet imaginary-time slices.
            // from the perspective of effiency, we donot need to calculate expV explicitly,
            // because in many cases, depending on the specific model and HS transformation,
            // they are diagonalized and it is straightforward to directly define their muplication rules.
            virtual void mult_B_from_left       ( GreensFunc& green, TimeIndex time_index, Spin spin ) = 0;
            virtual void mult_B_from_right      ( GreensFunc& green, TimeIndex time_index, Spin spin ) = 0;
            virtual void mult_invB_from_left    ( GreensFunc& green, TimeIndex time_index, Spin spin ) = 0;
            virtual void mult_invB_from_right   ( GreensFunc& green, TimeIndex time_index, Spin spin ) = 0;
            virtual void mult_transB_from_left  ( GreensFunc& green, TimeIndex time_index, Spin spin ) = 0;
    };

}

#endif // MODEL_BASE_H