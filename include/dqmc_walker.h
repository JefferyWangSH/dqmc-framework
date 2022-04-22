#ifndef DQMC_WALKER_H
#define DQMC_WALKER_H
#pragma once


/**
  *  This header file defines the crucial class QuantumMonteCarlo::DqmcWalker 
  *  to organize the entire dqmc program.
  */

#include <memory>
#include <vector>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace Utils { class SvdStack; }
namespace Model { class ModelBase; }
namespace Lattice { class LatticeBase; }
namespace Measure { class MeasureHandler; }


namespace QuantumMonteCarlo {

    // forward declaration
    class DqmcInitializer;

    using LatticeBase = Lattice::LatticeBase;
    using ModelBase = Model::ModelBase;
    using MeasureHandler = Measure::MeasureHandler;


    // ---------------------------- Crucial class QuantumMonteCarlo::DqmcWalker ----------------------------   
    class DqmcWalker {
        private:

            using TimeIndex = int;
            using RealScalar = double;
            using ptrRealScalarVec = std::unique_ptr<std::vector<double>>;
            using SvdStack = Utils::SvdStack;
            using ptrSvdStack = std::unique_ptr<SvdStack>;

            using GreensFunc = Eigen::MatrixXd;
            using GreensFuncVec = std::vector<Eigen::MatrixXd>;
            using ptrGreensFunc = std::unique_ptr<Eigen::MatrixXd>;
            using ptrGreensFuncVec = std::unique_ptr<std::vector<Eigen::MatrixXd>>;

            
            // --------------------------------- Walker params ---------------------------------------------

            int m_space_size{};                  // number of space sites
            int m_time_size{};                   // number of time slices
            RealScalar m_beta{};                 // inverse temperature beta
            RealScalar m_time_interval{};        // interval of imaginary-time grids
            int m_current_time_slice{};          // helping params to record current time slice


            // ----------------------- ( Equal-time and dynamic ) Greens functions -------------------------

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

            bool m_is_equaltime{};
            bool m_is_dynamic{};


            // ------------------------- SvdStack for numerical stabilization ------------------------------

            // Utils::SvdStack class
            // for efficient svd decompositions and numerical stabilization
            ptrSvdStack m_svd_stack_left_up{};
            ptrSvdStack m_svd_stack_left_dn{};
            ptrSvdStack m_svd_stack_right_up{};
            ptrSvdStack m_svd_stack_right_dn{};

            // pace of numerical stabilizations
            // or equivalently, the number of consequent wrapping steps of equal-time greens functions
            int m_stabilization_pace{};

            // keep track of the wrapping error
            RealScalar m_wrap_error{};


            // ---------------------------------- Reweighting params ---------------------------------------
            // keep track of the sign problem
            RealScalar m_config_sign{};
            ptrRealScalarVec m_vec_config_sign{};


        public:

            DqmcWalker() = default;

            // -------------------------------- Interfaces and friend class --------------------------------

            const int TimeSliceNum() const          { return this->m_time_size; }
            const RealScalar Beta()  const          { return this->m_beta; }
            const RealScalar TimeInterval() const   { return this->m_time_interval; }
            const RealScalar WrapError() const      { return this->m_wrap_error; }
            const int StabilizationPace() const     { return this->m_stabilization_pace; }

            // interface for greens functions
            // todo: this may cause problems if the pointer is nullptr
            GreensFunc& GreenttUp() { return *this->m_green_tt_up; }
            GreensFunc& GreenttDn() { return *this->m_green_tt_dn; }
            GreensFunc& Greent0Up() { return *this->m_green_t0_up; }
            GreensFunc& Greent0Dn() { return *this->m_green_t0_dn; }
            GreensFunc& Green0tUp() { return *this->m_green_0t_up; }
            GreensFunc& Green0tDn() { return *this->m_green_0t_dn; }

            GreensFuncVec& vecGreenttUp() { return *this->m_vec_green_tt_up; }
            GreensFuncVec& vecGreenttDn() { return *this->m_vec_green_tt_dn; }
            GreensFuncVec& vecGreent0Up() { return *this->m_vec_green_t0_up; }
            GreensFuncVec& vecGreent0Dn() { return *this->m_vec_green_t0_dn; }
            GreensFuncVec& vecGreen0tUp() { return *this->m_vec_green_0t_up; }
            GreensFuncVec& vecGreen0tDn() { return *this->m_vec_green_0t_dn; }

            friend class DqmcInitializer;
            

            // ------------------------------- Setup of parameters -----------------------------------------

            // set up the physical parameters
            // especially the inverse temperature and the number of time slices
            void set_physical_params( RealScalar beta, int time_size );

            // set up the pace of stabilizations
            void set_stabilization_pace( int stabilization_pace );


        private:

            // ---------------------------------- Initializations ------------------------------------------
            // never explicitly call these functions to avoid unpredictable mistakes,
            // and use DqmcInitializer instead 

            void initial( const LatticeBase& lattice, const MeasureHandler& meas_handler );

            void initial_svd_stacks( const LatticeBase& lattice, const ModelBase& model );

            // caution that this is a member function to initialize the model module
            // svd stacks should be initialized in advance
            void initial_greens_functions();

            // compute the sign of the initial bosonic configurations
            void initial_config_sign();

            // allocate memory
            void allocate_svd_stacks();
            void allocate_greens_functions();

        
        public:

            // ---------------------------------- Monte Carlo updates --------------------------------------
            
            // sweep forwards from time slice 0 to beta
            void sweep_from_0_to_beta( ModelBase& model );

            // sweep backwards from time slice beta to 0
            void sweep_from_beta_to_0( ModelBase& model );

            // sweep backwards from beta to 0 especially to compute dynamic greens functions,
            // without the updates of bosonic fields
            void sweep_for_dynamic_greens( ModelBase& model );

            
        private:

            // update the bosonic fields at time slice t using Metropolis algorithm
            void metropolis_update( ModelBase& model, TimeIndex t );
            
            // wrap the equal-time greens functions from time slice t to t+1
            void wrap_from_0_to_beta( const ModelBase& model, TimeIndex t );

            // wrap the equal-time greens functions from time slice t to t-1
            void wrap_from_beta_to_0( const ModelBase& model, TimeIndex t );

    };

}


#endif // DQMC_WALKER_H