#ifndef DQMC_HUBBARD_DETQMC_H
#define DQMC_HUBBARD_DETQMC_H
#pragma once

/**
 *  This head file includes detQMC class
 *  for the Monte Carlo simulation of half-filled Hubbard model.
 *  including: MC sampling (measuring) and output of statistical results.
 */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2

#include <map>
#include "hubbard.h"


class detQMC
{
private:

    // model parameters
    Hubbard hubb;
    int nwrap{10}, nwarm{300};
    int nbin{20}, nsweep{100}, nBetweenBins{10};

    // bool parameters to control thermalization and measurements
    bool bool_warm_up = true;
    bool bool_measure_eqtime = true;
    bool bool_measure_dynamic = true;

    // for equal-time measurements
    std::map<std::string, std::vector<double>> obs_bin_eqtime;
    std::map<std::string, double> obs_mean_eqtime;
    std::map<std::string, double> obs_err_eqtime;

    // for time-displaced measurements: Matsubara Green's function gt0
    matXd obs_bin_gt0[500][500];        // obs_bin_gt0 [bin][tau] <type matXd>
    double obs_mean_gt0_k[500];         // mean value: obs_mean_gt0_k[tau] <double>
    double obs_err_gt0_k[500];          // statistical error: obs_err_gt0_k[tau] <double>

    /* temporary parameters for measurements */
    // equal-time
    int n_equal_time = 0;
    double DoubleOccu = 0.0;
    double KineticEnergy = 0.0;
    double StructFactor = 0.0;
    double MomentumDist = 0.0;
    double localSpinCorr = 0.0;

    // time-displaced
    int n_time_displaced = 0;
    matXd vec_gt0_tau[500];

    // lattice momentum q
    vecXd q = vecXd::Zero(2);

    // time cost of one single measuring process
    time_t  begin_t, end_t;


public:

    detQMC() = default;

    ~detQMC();

    /* set up model parameters */
    void set_Model_Params(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap);

    /* set up parameters for Monte Carlo simulation */
    void set_MC_Params(int nwarm, int nbin, int nsweep, int nBetweenBins);

    /* set up bool parameters */
    void set_bool_Params(bool bool_warm_up, bool bool_measure_eqtime, bool bool_measure_dynamic);

    /* FIXME: set up momentum q, modified to allow sequence of q */
    void set_Momentum_q(double qx, double qy);

    /* prepare for measuring */
    void initialMeasure();

    /* run a dqmc simulation */
    void runQMC(bool bool_display_process);

    /* analyse statistics from simulation */
    void analyseStats();

    /* print results of measurements in the command console */
    void printStats();

    /* write results of measurements to file */
    void output_Stats_eqtime(const std::string& filename, bool bool_Append);

    void output_Stats_dynamic(const std::string& filename, bool bool_Append);

private:

    /* process of back-and-forth sweep, meantime do the measurements */
    void sweep_BackAndForth(bool bool_eqtime, bool bool_dynamic);

    /* equal-time measurements */
    void measure_equal_time();

    /* time-displaced measurements */
    void measure_time_displaced();

    /** double occupation: D = < n_up*n_dn > */
    void meas_Double_Occu(const matXd& gu, const matXd& gd);

    /** single particle kinetic energy */
    void meas_Kinetic_Energy(const matXd& gu, const matXd& gd);

    /** momentum distribution of electrons: fourier transformation of real-space electron distribution */
    void meas_Momentum_Dist(const matXd& gu, const matXd& gd, const vecXd& p);

    /** local spin correlation: magnetization C(0,0) = < (n_up - n_dn)^2 > */
    void meas_local_Spin_Corr(const matXd& gu, const matXd& gd);

    /** magnetic struct factor: fourier transformation of real-space spin-spin correlation */
    void meas_Struct_Factor(const matXd& gu, const matXd& gd, const vecXd& p);

    /* normalize data from scratch */
    void normalizeStats();

    /* bin measurements */
    void write_Stats_to_bins(int bin);

    /* analyse equal-time data and compute means and errors */
    void analyse_equal_time_Stats(const std::string& obs);

    /** time-displaced green function in momentum space < c(k,tau) * c^T(k,0) > */
    void analyse_timeDisplaced_Stats(const vecXd& p);

    /* clear data in memory */
    void clearStats();
};

#endif //DQMC_HUBBARD_DETQMC_H
