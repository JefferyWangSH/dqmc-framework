#ifndef DQMC_HUBBARD_MEASURE_METHODS_H
#define DQMC_HUBBARD_MEASURE_METHODS_H
#pragma once

/**
  *  This header file includes declarations of user-defined methods 
  *  for the measurements of physical observables using QMC.
  *  Both equal-time and dynamical measurements are supported.
  */

#include "observable.h"

// forward declaration
namespace Model { class Hubbard; }

namespace Measure {

    class Methods {
        public:
        // definitions of measuring methods
        // arguments of method functions should include Observable, Measure and Hubbard class

        // Equal-time Measurements:
        //    1. Filling number <n>
        //    2. Double occupancy D = < n_up*n_dn >
        //    3. Single particle kinetic energy
        //    4. Distributions of electrons in momentum space
        //    5. Local spin correlation, magnetization C(0,0) = < (n_up - n_dn)^2 >
        //    6. Spin density structure factor (SDW)
        //    7. Charge density structure factor (CDW)
        //    8. S wave Cooper pairing correlation function

        static void measure_config_sign_eqtime(Observable<double> &sign_eqtime, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_filling_number(Observable<double> &filling_number, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_double_occupancy(Observable<double> &double_occupancy, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_kinetic_energy(Observable<double> &kinetic_energy, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_local_spin_corr(Observable<double> &local_spin_corr, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_momentum_distribution(Observable<double> &momentum_dist, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_spin_density_structure_factor(Observable<double> &sdw_factor, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_charge_density_structure_factor(Observable<double> &cdw_factor, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_s_wave_pairing_corr(Observable<double> &s_wave_pairing, Measure &measure, const Model::Hubbard &hubbard);


        // Dynamical Measurements:
        //    1. Dynamical green's functions in momentum space: G(k, tau) = < c(k, tau) * c^+(k, 0) >
        //    2. Superfluid stiffness \rho_s of superconducting: \rho_s = (\Gamma_L - \Gamma_T) / 4
        //    3. Density of states in imaginary time space: N(tau) = 1/N * \sum_{i} G(tau, 0)_{ii}
        
        static void measure_config_sign_dynamic(Observable<double> &sign_eqtime, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_greens_functions(Observable<Eigen::MatrixXd> &greens_functions, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_density_of_states(Observable<Eigen::VectorXd> &density_of_states, Measure &measure, const Model::Hubbard &hubbard);
        static void measure_superfluid_stiffness(Observable<double> &density_of_states, Measure &measure, const Model::Hubbard &hubbard);
    };

} // namespace Measure

#endif //DQMC_HUBBARD_MEASURE_METHODS_H