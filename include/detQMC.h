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
#include "hubbard.h"

#include <map>
#include <chrono>


class detQMC {

private:
    // model parameters
    Hubbard hubb;
    int nwrap{10}, nwarm{300};
    int nbin{20}, nsweep{100}, nBetweenBins{10};

    // bool parameters to control thermalization and measurements
    bool bool_warm_up = true;
    bool bool_measure_eqtime = true;
    bool bool_measure_dynamic = true;

    // lattice momentum q
    Eigen::VectorXd q = Eigen::VectorXd::Zero(2);

    // time cost of one single measuring process
    std::chrono::steady_clock::time_point begin_t{}, end_t{};


public:

    // for equal-time measurements
    measure::eqtimeMeasure eqtimeMeasure;

    // for time-displaced (dynamic) measurements
    measure::dynamicMeasure dynamicMeasure;

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

    /* print out simulation params on the command console */
    void printParams();

    /* prepare for measuring */
    void initialMeasure();

    /* run a dqmc simulation */
    void runQMC(bool bool_display_process);

    /* analyse statistics from simulation */
    void analyseStats();

    /* print results of measurements on the command console */
    void printStats();

    /* write results of measurements to file */
    void output_Stats_eqtime(const std::string& filename, bool bool_Append);

    void output_Stats_dynamic(const std::string& filename, bool bool_Append);

private:

    /* process of back-and-forth sweep, meantime do the measurements */
    void sweep_BackAndForth(bool bool_eqtime, bool bool_dynamic);
};

#endif //DQMC_HUBBARD_DETQMC_H
