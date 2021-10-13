#ifndef DQMC_HUBBARD_DETQMC_H
#define DQMC_HUBBARD_DETQMC_H
#pragma once

/**
  *  This head file includes detQMC class
  *  for the Monte Carlo simulation of half-filled Hubbard model.
  *  including: MC sampling (measuring) and output of statistical results.
  */

namespace Model { class Hubbard; }
namespace Measure { class EqtimeMeasure; }
namespace Measure { class DynamicMeasure; }
class SvdStack;

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <chrono>


namespace Simulation {

    class DetQMC {
    private:
        // model parameters
        Model::Hubbard *hubb{};

        int nwrap{10}, nwarm{300};
        int nbin{20}, nsweep{100}, n_between_bins{10};

        // bool parameters to control thermalization and measurements
        bool bool_warm_up = true;
        bool bool_measure_eqtime = true;
        bool bool_measure_dynamic = true;

        // lattice momentum q
        Eigen::VectorXd q = Eigen::VectorXd::Zero(2);

        // time cost of one single measuring process
        std::chrono::steady_clock::time_point begin_t{}, end_t{};

        // for equal-time measurements
        Measure::EqtimeMeasure *EqtimeMeasure{};

        // for time-displaced (dynamical) measurements
        Measure::DynamicMeasure *DynamicMeasure{};

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        DetQMC() = default;

        ~DetQMC();

        /** Following functions aimed to set up params and initialization before DQMC calculation starts */
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


        /** Critical functions for DQMC calculations, including Monte Carlo updates and measurements **/
        /* initialization, especially preparing for measuring */
        void init_measure();

        /* run a dqmc simulation */
        void run_QMC(bool bool_display_process);

        /* analyse statistics from simulation */
        void analyse_stats() const;


        /** Output modules for the output of DQMC measuring results */
        /* print simulation params onto terminal */
        void print_params() const;

        /* print measuring results onto terminal */
        void print_stats() const;

        /* write results of measurements, including means and errors, into file */
        void file_output_eqtime_stats(const std::string &filename) const;

        void file_output_dynamic_stats(const std::string &filename) const;

        /* write sequences of imaginary-time tau into file */
        void file_output_tau(const std::string &filename) const;

        /* write results of green's function measurements, in terms of bins, into file */
        void bin_output_greens(const std::string &filename) const;

        /* write measuring results of local density of states(LDOS) , in terms of bins, into file */
        void bin_output_LDOS(const std::string &filename) const;

        /* analyse space correlation of cooper pairs and write results into file */
        void file_output_cooper_corr(const std::string &filename) const;

        /* output aux field configurations into file */
        void file_output_aux_field_configs(const std::string &filename) const;

    private:
        /* process of Monte Carlo sweeping, and do the measurements if needed */
        void sweep_back_and_forth(bool bool_eqtime, bool bool_dynamic) const;
    };

}


#endif //DQMC_HUBBARD_DETQMC_H
