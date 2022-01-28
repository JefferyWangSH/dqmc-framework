#ifndef DQMC_HUBBARD_DETQMC_H
#define DQMC_HUBBARD_DETQMC_H
#pragma once

/**
  *  This head file includes detQMC class
  *  for the Monte Carlo simulation of half-filled Hubbard model.
  *  including: MC sampling (measuring) and output of statistical results.
  */

#include <chrono>
#include <memory>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "hubbard.h"
#include "measure.h"

namespace FileOutput { 
    void file_output_tau(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode); 
    void file_output_aux_field(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode); 
}
namespace ScreenOutput {
    void screen_output_params(const int &world_size, const Simulation::DetQMC &dqmc); 
    void screen_output_end_info(const Simulation::DetQMC &dqmc);
}


namespace Simulation {

    class DetQMC {
    private:
        // model parameters
        std::unique_ptr<Model::Hubbard> hubbard{};

        int nwrap{10}, nwarm{300};
        int nbin{20}, nsweep{100}, n_between_bins{10};

        // bool parameters to control thermalization and measurements
        bool is_warm_up{true};
        bool is_eqtime_measure{true};
        bool is_dynamic_measure{true};

        // input configuration file and list of observables
        std::unique_ptr<std::string> config_file{};
        std::unique_ptr<std::vector<std::string>> obs_list{};

        // lattice momentum q
        Eigen::VectorXd q = Eigen::VectorXd::Zero(2);

        // time cost of one single measuring process
        std::chrono::steady_clock::time_point begin_t{}, end_t{};

    public:
        // QMC measurements
        std::unique_ptr<Measure::Measure> measure{};

        // friend function
        friend void FileOutput::file_output_tau(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode);
        friend void FileOutput::file_output_aux_field(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode); 
        friend void ScreenOutput::screen_output_params(const int &world_size, const Simulation::DetQMC &dqmc);
        friend void ScreenOutput::screen_output_end_info(const Simulation::DetQMC &dqmc);

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        DetQMC() = default;
        ~DetQMC() = default;

        /** functions to set up params and initialization before DQMC calculation starts */
        /* set up model parameters */
        void set_model_params(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap);

        /* set up parameters for Monte Carlo simulation */
        void set_Monte_Carlo_params(int nwarm, int nbin, int nsweep, int n_between_bins);

        /* set up bool parameters */
        void set_controlling_params(bool is_warm_up, bool is_eqtime_measure, bool is_dynamic_measure, bool is_checkerboard);

        /* set up list of observables to be measured */
        void set_observable_list(const std::vector<std::string> &obs_list);

        /* set up input file of aux field configurations */
        void set_aux_field_configs(const std::string &config_file);

        /* set up lattice momentum q for measurements in momentum space */
        void set_lattice_momentum(double qx, double qy);

        /** Critical functions for DQMC calculations, including Monte Carlo updates and measurements **/
        /* initialization, especially allocating memory for monte carlo and measurements */
        void initial();

        /* run a dqmc simulation */
        void run(bool show_running_process);

        /* analyse statistics from simulation */
        void analyse_stats() const;

    private:
        /* process of Monte Carlo sweeping, and do the measurements if needed */
        void sweep_forth_and_back(bool is_eqtime_measure, bool is_dynamic_measure) const;

        /* read input configurations from input file */
        void read_configs_from_file(const std::string &config_file);

        /* initialize hubbard class with input configs */
        void initial_hubbard_with_input_configs();
    };

}

#endif //DQMC_HUBBARD_DETQMC_H