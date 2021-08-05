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
    int nbin{20}, nsweep{100}, n_between_bins{10};
    bool is_checkerboard = true;

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
    void set_model_params(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap, bool is_checkerboard);

    /* set up parameters for Monte Carlo simulation */
    void set_Monte_Carlo_params(int nwarm, int nbin, int nsweep, int n_between_bins);

    /* set up bool parameters */
    void set_controlling_params(bool bool_warm_up, bool bool_measure_eqtime, bool bool_measure_dynamic);

    /* set up lattice momentum q for momentum measurements */
    void set_lattice_momentum(double qx, double qy);

    /* read aux field configurations from input file */
    void read_aux_field_configs(const std::string &filename);

    /* print out simulation params on the command console */
    void print_params();

    /* prepare for measuring */
    void init_measure();

    /* run a dqmc simulation */
    void run_QMC(bool bool_display_process);

    /* analyse statistics from simulation */
    void analyse_stats();

    /* print results of measurements on the command console */
    void print_stats();

    /* write results of measurements to file */
    void file_output_stats_eqtime(const std::string &filename, bool bool_Append);

    void file_output_stats_dynamic(const std::string &filename, bool bool_Append);

    /* output aux field configurations to file */
    void file_output_aux_field_configs(const std::string &filename);

public:

    /* process of back-and-forth sweep, meantime do the measurements */
    void sweep_back_and_forth(bool bool_eqtime, bool bool_dynamic);
};

#endif //DQMC_HUBBARD_DETQMC_H
