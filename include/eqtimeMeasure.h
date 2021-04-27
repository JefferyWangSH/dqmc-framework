#ifndef DQMC_HUBBARD_EQTIMEMEASURE_H
#define DQMC_HUBBARD_EQTIMEMEASURE_H
#pragma once

/**
 *  This head file includes module for equal-time measuring.
 *  Class: measure::eqtimeMeasure
 *  Measuring:
 *   1. Double occupancy D = < n_up*n_dn >
 *   2. Single particle kinetic energy
 *   3. Density-density correlation in momentum space
 *   4. Local spin correlation, magnetization C(0,0) = < (n_up - n_dn)^2 >
 *   5. Magnetic struct factor, spin-spin correlation in momentum space
 *   6. ...
 */

#include <map>
#include <vector>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

// forward declaration
class Hubbard;

namespace measure{
    class eqtimeMeasure {
    public:
        int nbin{20};

        // temporary parameters
        int n_equal_time = 0;
        double DoubleOccu = 0.0;
        double KineticEnergy = 0.0;
        double StructFactor = 0.0;
        double MomentumDist = 0.0;
        double localSpinCorr = 0.0;

        // for equal-time measurements
        std::map<std::string, std::vector<double>> obs_bin_eqtime;
        std::map<std::string, double> obs_mean_eqtime;
        std::map<std::string, double> obs_err_eqtime;

        // lattice momentum q
        Eigen::VectorXd q = Eigen::VectorXd::Zero(2);


        /* construct function */
        eqtimeMeasure()  = default;

        /* resize the size of bins */
        void resize(const int &nbin);

        /* prepare for measuring */
        void initial();

        /* clear temporary parameters */
        void clear();

        /* equal-time measurements */
        void measure_equal_time(const Hubbard &hubbard);

        /* normalize data from scratch */
        void normalizeStats(const Hubbard &hubbard);

        /* bin measurements */
        void write_Stats_to_bins(int bin);

        /* analyse certain equal-time data and compute means and errors */
        void analyse_equal_time_Stats(const std::string& obs);

        /* analyse all equal-time statistics from simulation */
        void analyseStats();

    private:
        /** double occupation: D = < n_up*n_dn > */
        void meas_Double_Occu(const Hubbard &hubbard, const int &t);

        /** single particle kinetic energy */
        void meas_Kinetic_Energy(const Hubbard &hubbard, const int &t);

        /** momentum distribution of electrons: fourier transformation of real-space electron distribution */
        void meas_Momentum_Dist(const Hubbard &hubbard, const int &t, const Eigen::VectorXd& p);

        /** local spin correlation: magnetization C(0,0) = < (n_up - n_dn)^2 > */
        void meas_local_Spin_Corr(const Hubbard &hubbard, const int &t);

        /** magnetic struct factor: fourier transformation of real-space spin-spin correlation */
        void meas_Struct_Factor(const Hubbard &hubbard, const int &t, const Eigen::VectorXd& p);
    };
}

#endif //DQMC_HUBBARD_EQTIMEMEASURE_H
