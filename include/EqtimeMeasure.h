#ifndef DQMC_HUBBARD_EQTIMEMEASURE_H
#define DQMC_HUBBARD_EQTIMEMEASURE_H
#pragma once

/**
  *  This head file includes module for equal-time measuring.
  *  Class: measure::eqtimeMeasure
  *  Measuring:
  *   1. Double occupancy D = < n_up*n_dn >
  *   2. Single particle kinetic energy
  *   3. Density of electrons in momentum space
  *   4. Local spin correlation, magnetization C(0,0) = < (n_up - n_dn)^2 >
  *   5. Magnetic struct factor, that is, spin-spin correlation in momentum space
  *   6. Space correlation of Copper-like order parameter
  *   7. ...
  */

#include <map>
#include <vector>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "MeasureData.h"

// forward declaration
namespace Model { class Hubbard; }
namespace Measure { class MeasureData; }


namespace Measure{

    class EqtimeMeasure {
    public:
        int nbin{20};

        /* for equal-time (static) measurements */

        Measure::MeasureData double_occu;                 // double occupancy
        Measure::MeasureData kinetic_energy;              // kinetic energy
        Measure::MeasureData electron_density;            // electron density in momentum space
        Measure::MeasureData local_corr;                  // local spin correlation
        Measure::MeasureData AFM_factor;                  // AFM structure factor
        std::vector<Measure::MeasureData> cooper_corr;    // space correlation of Cooper-type order parameter
        Measure::MeasureData sign;                        // average sign to keep track of sign problem

        // measurements of equal-time greens functions
        std::vector<Eigen::MatrixXd> bin_gtt_up;
        std::vector<Eigen::MatrixXd> bin_gtt_dn;

        // temporary parameters
        int n_equal_time = 0;
        double tmp_sign = 0.0;
        Eigen::MatrixXd tmp_gtt_up;
        Eigen::MatrixXd tmp_gtt_dn;

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
        void measure_equal_time_greens(const Model::Hubbard &hubbard);

        /* normalize data from scratch */
        void normalizeStats(const Model::Hubbard &hubbard);

        /* bin measurements */
        void write_Stats_to_bins(int bin);

//        /* analyse certain equal-time data and compute means and errors */
//        void analyse_equal_time_Stats(const std::string &obs);

        /* analyse all equal-time statistics from simulation */
        void analyseStats(const Model::Hubbard &hubbard);

    private:
        /** double occupation: D = < n_up*n_dn > */
        void analyse_double_occu(const int &bin, const Model::Hubbard &hubbard);

        /** single particle kinetic energy */
        void analyse_kinetic_energy(const int &bin, const Model::Hubbard &hubbard);

        /** momentum distribution of electrons: fourier transformation of real-space electron distribution */
        void analyse_electron_density(const int &bin, const Model::Hubbard &hubbard);

        /** local spin correlation: magnetization C(0,0) = < (n_up - n_dn)^2 > */
        void analyse_local_corr(const int &bin, const Model::Hubbard &hubbard);

        /** Anti-ferromagnetic structure factor: fourier transformation of real-space pi-pi correlation of spins */
        void analyse_AFM_factor(const int &bin, const Model::Hubbard &hubbard);

        /** space correlation of Cooper-type order parameter */
        void analyse_Cooper_corr(const int &bin, const Model::Hubbard &hubbard);
    };
}

#endif //DQMC_HUBBARD_EQTIMEMEASURE_H
