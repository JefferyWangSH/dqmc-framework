#ifndef DQMC_HUBBARD_DYNAMICMEASURE_H
#define DQMC_HUBBARD_DYNAMICMEASURE_H
#pragma once

/**
  *  This head file includes module for time-displaced (dynamical) measuring.
  *  Class: Measure::dynamicMeasure
  *  Measuring:
  *   1. Dynamical green's function of imaginary time: G(k, tau) = < c(k, tau) * c^+(k, 0) >
  *   2. Superfluid stiffness \rho_s of superconducting: \rho_s = (\Gamma_L - \Gamma_T) / 4
  *   3. Density of states in imaginary time space: N(tau) = 1/N * \sum_{i} G(tau, 0)_{ii}
  *   4. ...
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <vector>
#include "observable.h"


// forward declaration
namespace Model { class Hubbard; }

namespace Measure{
    class DynamicMeasure {
    public:
        int nbin{20};

        /* for time-displaced measurements */
        // dynamical correlation function of imaginary time G(k, \tua) = < c(k, \tau) * c^+(k, 0) >, with \tau > 0.
        Measure::Observable<Eigen::VectorXd> matsubara_greens;

        // density of states (DOS) measurements 1/N * \sum_{i} G(\tau, 0)_{ii} or G(\tau, 0).trace()/N
        Measure::Observable<Eigen::VectorXd> density_of_states;

        // superfluid stiffness (helicity modulus) \rho_s of superconducting
        Measure::Observable<double> superfluid_stiffness;
        // Measure::Observable<Eigen::MatrixXd> current_current_corr;

        // sign problem
        Measure::Observable<double> sign;

        // lattice momentum q
        Eigen::VectorXd q = Eigen::VectorXd::Zero(2);

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /* (de)construct functions */
        DynamicMeasure() = default;
        explicit DynamicMeasure(const int &nbin);

        ~DynamicMeasure() = default;

        /* resize the size of bins */
        void resize(const int &nbin);

        /* prepare for measuring */
        void initial(const Model::Hubbard &hubbard);

        /* clear temporary parameters */
        void clear_temporary(const Model::Hubbard &hubbard);

        /* time-displaced measurements */
        void time_displaced_measure(const Model::Hubbard &hubbard);

        /* normalize data from scratch */
        void normalize_stats(const Model::Hubbard &hubbard);

        /* bin measurements */
        void write_stats_to_bins(const int &bin, const Model::Hubbard &hubbard);

        /* analyse dynamical statistics */
        void analyse_stats(const Model::Hubbard &hubbard);

    private:
        /** dynamical correlation function in momentum space < c(k,tau) * c^+(k,0) > */
        void measure_matsubara_greens(const int &t, const Model::Hubbard &hubbard);

        /** density of states (in imaginary time) 1/N * \sum_{i} G(tau, 0)_{ii} */
        void measure_density_of_states(const int &t, const Model::Hubbard &hubbard);

        /** superfluid stiffness (helicity modules) \rho_s */
        void measure_superfluid_stiffness(const Model::Hubbard &hubbard);
    };
}

#endif //DQMC_HUBBARD_DYNAMICMEASURE_H
