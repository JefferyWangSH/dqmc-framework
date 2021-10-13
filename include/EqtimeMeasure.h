#ifndef DQMC_HUBBARD_EQTIMEMEASURE_H
#define DQMC_HUBBARD_EQTIMEMEASURE_H
#pragma once

/**
  *  This head file includes module for equal-time measuring.
  *  Class: Measure::eqtimeMeasure
  *  Measuring:
  *   1. Double occupancy D = < n_up*n_dn >
  *   2. Single particle kinetic energy
  *   3. Density of electrons in momentum space
  *   4. Local spin correlation, magnetization C(0,0) = < (n_up - n_dn)^2 >
  *   5. Magnetic struct factor, that is, spin-spin correlation in momentum space
  *   6. Space correlation of Copper-like order parameter
  *   7. ...
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <vector>
#include "MeasureData.h"

// forward declaration
namespace Model { class Hubbard; }
namespace Measure { class MeasureData; }


namespace Measure{

    class EqtimeMeasure {
    public:
        int nbin{20};

        /* for equal-time (static) measurements */
        Measure::MeasureData sign;                                      // average sign to keep track of sign problem
        Measure::MeasureData double_occu;                               // double occupancy
        Measure::MeasureData kinetic_energy;                            // kinetic energy
        Measure::MeasureData electron_density;                          // electron density in momentum space
        Measure::MeasureData local_corr;                                // local spin correlation
        Measure::MeasureData AFM_factor;                                // AFM structure factor
        Eigen::MatrixX<Measure::MeasureData> cooper_corr;               // space correlation of (local) Cooper order parameter

        // temporary counting parameters
        int n_equal_time = 0;

        // lattice momentum q
        Eigen::VectorXd q = Eigen::VectorXd::Zero(2);

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /* (de)construct functions */
        EqtimeMeasure()  = default;
        explicit EqtimeMeasure(const int &nbin);

        ~EqtimeMeasure() = default;

        /* resize the size of bins */
        void resize(const int &nbin);

        /* prepare for measuring */
        void initial(const Model::Hubbard &hubbard);

        /* clear temporary parameters */
        void clear_temporary();

        /* equal-time measurements */
        void equal_time_measure(const Model::Hubbard &hubbard);

        /* normalize data from scratch */
        void normalize_stats(const Model::Hubbard &hubbard);

        /* bin measurements */
        void write_stats_to_bins(int bin);

        /* analyse all equal-time statistics from simulation */
        void analyse_stats(const Model::Hubbard &hubbard);

    private:
        /** double occupation: D = < n_up * n_dn > */
        void measure_double_occu(const int &tau, const Model::Hubbard &hubbard);

        /** single particle kinetic energy */
        void measure_kinetic_energy(const int &tau, const Model::Hubbard &hubbard);

        /** momentum distribution of electrons: fourier transformation of real-space electron distribution */
        void measure_electron_density(const int &tau, const Model::Hubbard &hubbard);

        /** local spin correlation: magnetization C(0,0) = < (n_up - n_dn)^2 > */
        void measure_local_corr(const int &tau, const Model::Hubbard &hubbard);

        /** Anti-ferromagnetic structure factor: fourier transformation of real-space pi-pi correlation of spins */
        void measure_AFM_factor(const int &tau, const Model::Hubbard &hubbard);

        /** space correlation of (local) Cooper order parameter */
        void measure_Cooper_corr(const int &tau, const Model::Hubbard &hubbard);
    };
}

#endif //DQMC_HUBBARD_EQTIMEMEASURE_H
