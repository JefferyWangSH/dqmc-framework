#ifndef DQMC_WALKER_H
#define DQMC_WALKER_H
#pragma once


/**
  *  This header file defines the crucial class QuantumMonteCarlo::DqmcWalker 
  *  to organize the entire dqmc program.
  */

#include <memory>
#include <vector>

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


    // ------------------------ Crucial class QuantumMonteCarlo::DqmcWalker ----------------------------   
    class DqmcWalker {
        private:

            using TimeIndex = int;
            using RealScalar = double;
            using ptrRealScalarVec = std::unique_ptr<std::vector<double>>;
            using SvdStack = Utils::SvdStack;
            using ptrSvdStack = std::unique_ptr<SvdStack>;
            
            int m_space_size{};                  // number of space sites
            int m_time_size{};                   // number of time slices
            RealScalar m_beta{};                 // inverse temperature beta
            RealScalar m_time_interval{};        // interval of imaginary-time grids
            int m_current_time_slice{};          // helping params to record current time slice

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

            // keep track of the sign problem
            RealScalar m_config_sign{};
            ptrRealScalarVec m_vec_config_sign{};

            bool m_is_equaltime{};
            bool m_is_dynamic{};

        public:

            DqmcWalker() = default;

            // --------------------------- Interfaces and friend class -------------------------

            const int TimeSliceNum() const;
            const RealScalar Beta()  const;
            const RealScalar TimeInterval() const;
            const RealScalar WrapError() const;
            const int StabilizationPace() const;

            friend class DqmcInitializer;
            

            // ------------------------------ Setup of parameters ------------------------------

            // set up the physical parameters
            // especially the inverse temperature and the number of time slices
            void set_physical_params( RealScalar beta, int time_size );

            // set up the pace of stabilizations
            void set_stabilization_pace( int stabilization_pace );


        private:
        
            // ------------------------------- Initialization ----------------------------------
            // never explicitly call these functions to avoid unpredictable mistakes, and use DqmcInitializer instead 

            void initial( const LatticeBase& lattice, const MeasureHandler& meas_handler );

            void initial_svd_stacks( const LatticeBase& lattice, const ModelBase& model );

            // caution that this is a member function to initialize the model module
            // svd stacks should be initialized in advance
            void initial_greens_function( ModelBase& model );

        
        public:

            // ----------------------------- Monte Carlo updates -------------------------------
            
            // sweep forwards from time slice 0 to beta
            void sweep_from_0_to_beta( ModelBase& model );

            // sweep backwards from time slice beta to 0
            void sweep_from_beta_to_0( ModelBase& model );

            // sweep backwards from beta to 0 especially to compute dynamic greens functions,
            // without the updates of bosonic fields
            void sweep_for_dynamic_greens( ModelBase& model );

            // wrap the equal-time greens functions from time slice t to t+1
            void wrap_from_0_to_beta( ModelBase& model, TimeIndex t );

            // wrap the equal-time greens functions from time slice t to t-1
            void wrap_from_beta_to_0( ModelBase& model, TimeIndex t );

    };

}


#endif // DQMC_WALKER_H