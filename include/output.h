#ifndef DQMC_HUBBARD_OUTPUT_H
#define DQMC_HUBBARD_OUTPUT_H
#pragma once

/**
  *  This head file includes declarations of subroutines for outputing DQMC measuring results,
  *  in form of outfile stream or terminal output.
  */

#include "observable.h"

// forward declaration
namespace Simulation { class DetQMC; }

namespace FileOutput{

    void file_output_observable(const Measure::Observable<double> &obs, const std::string &file_name, const int &mode);
    void file_output_observable(const Measure::Observable<Eigen::VectorXd> &obs, const std::string &file_name, const int &mode);
    void file_output_observable(const Measure::Observable<Eigen::MatrixXd> &obs, const std::string &file_name, const int &mode);

    void file_output_observable_bin(const Measure::Observable<double> &obs, const std::string &file_name, const int &mode);
    void file_output_observable_bin(const Measure::Observable<Eigen::VectorXd> &obs, const std::string &file_name, const int &mode);
    void file_output_observable_bin(const Measure::Observable<Eigen::MatrixXd> &obs, const std::string &file_name, const int &mode);

    void file_output_tau(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode);
    void file_output_qlist(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode);
    void file_output_aux_field(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode);

} // namespace FileOutput


namespace ScreenOutput {

    void screen_output_init_info(const std::string &master_proc_name, const int &world_size, const Simulation::DetQMC &dqmc);

    void screen_output_end_info(const Simulation::DetQMC &dqmc);

    void screen_output_time();

    void screen_output_mpi(const std::string &master_proc_name, const int &world_size);

    void screen_output_params(const int &world_size, const Simulation::DetQMC &dqmc);

    // only print out observables of double type
    void screen_output_observable(const Measure::Observable<double> &obs, const std::string &obs_name);

} // namespace ScreenOutput

#endif // DQMC_HUBBARD_OUTPUT_H