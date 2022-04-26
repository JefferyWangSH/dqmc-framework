#ifndef MEASURE_METHODS_H
#define MEASURE_METHODS_H
#pragma once

/**
  *  This header file includes declarations of user-defined methods 
  *  for the measurements of physical observables using dqmc.
  *  Both equal-time and dynamical measurements are supported.
  */

#include "measure/observable.h"

// forward declaration
namespace Model { class ModelBase; }
namespace Lattice { class LatticeBase; }
namespace QuantumMonteCarlo { class DqmcWalker; }

namespace Measure {

    using ModelBase = Model::ModelBase;
    using LatticeBase = Lattice::LatticeBase;
    using DqmcWalker = QuantumMonteCarlo::DqmcWalker;
    using ScalarObs = Observable::ScalarObs;
    using VectorObs = Observable::VectorObs;
    using MatrixObs = Observable::MatrixObs;


    // ------------------------------ Interface class Measure::Method ---------------------------------
    // provide user-defined measuring methods
    class Methods {


        public:
            // definitions of measuring methods
            // arguments of method functions should include Observable<ObsType>, 
            // Measure::MeasureHandler, Model::ModelBase class and Lattice::LatticeBase

            // Equal-time Measurements:
            //    1. Filling number <n>
            //    2. Double occupancy D = < n_up*n_dn >
            //    3. Single particle kinetic energy
            //    4. Distributions of electrons in momentum space
            //    5. Local spin correlation, magnetization C(0,0) = < (n_up - n_dn)^2 >
            //    6. Spin density structure factor (SDW)
            //    7. Charge density structure factor (CDW)
            //    8. S wave Cooper pairing correlation function

            static void measure_equaltime_config_sign           (  ScalarObs& equaltime_sign, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice );

            static void measure_filling_number                  (  ScalarObs& filling_number, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice );

            static void measure_double_occupancy                (  ScalarObs& double_occupancy, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model, 
                                                                   const LatticeBase& lattice );
                                                                  
            static void measure_kinetic_energy                  (  ScalarObs& kinetic_energy, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice );

            static void measure_local_spin_corr                 (  ScalarObs& local_spin_corr, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice ) {};

            static void measure_momentum_distribution           (  ScalarObs& momentum_dist, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice ) {};

            static void measure_spin_density_structure_factor   (  ScalarObs& sdw_factor, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice ) {};

            static void measure_charge_density_structure_factor (  ScalarObs& cdw_factor, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice ) {};

            static void measure_s_wave_pairing_corr             (  ScalarObs& s_wave_pairing,
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice ) {};


            // Dynamical Measurements:
            //    1. Dynamical green's functions in momentum space: G(k, tau) = < c(k, tau) * c^+(k, 0) >
            //    2. Superfluid stiffness \rho_s of superconducting: \rho_s = (\Gamma_L - \Gamma_T) / 4
            //    3. Density of states in imaginary time space: N(tau) = 1/N * \sum_{i} G(tau, 0)_{ii}
            
            static void measure_dynamic_config_sign             (  ScalarObs& dynamic_sign, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice );

            static void measure_greens_functions                (  MatrixObs& greens_functions, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice ) {};

            static void measure_density_of_states               (  VectorObs& density_of_states, 
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice ) {};

            static void measure_superfluid_stiffness            (  ScalarObs& density_of_states,
                                                                   const DqmcWalker& walker,
                                                                   const ModelBase& model,
                                                                   const LatticeBase& lattice ) {};
    
    };

} // namespace Measure

#endif // MEASURE_METHODS_H