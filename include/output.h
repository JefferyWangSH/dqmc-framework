#ifndef DQMC_HUBBARD_OUTPUT_H
#define DQMC_HUBBARD_OUTPUT_H
#pragma once

/**
  *  This head file includes declarations of subroutines for the output of DQMC measuring results,
  *  in form of outfile stream or terminal output.
  * 
  */

#include <vector>
#include <string>

namespace Measure { class MeasureData; }
namespace Simulation { class DetQMC; }


namespace FileOutput{

    void file_output_observable(const Measure::MeasureData &obs, const std::string &file_name, const int &mode);

    void file_output_observable_bin(const Measure::MeasureData &obs, const std::string &file_name, const int &mode);

    void file_output_observable(const std::vector<Measure::MeasureData> &obs_vec, const std::string &file_name, const int &mode);

    void file_output_observable_bin(const std::vector<Measure::MeasureData> &obs_vec, const std::string &file_name, const int &mode);

} // namespace FileOutput


namespace ScreenOutput {

    void screen_output_init_info(const std::string &master_proc_name, const int &world_size, const Simulation::DetQMC &dqmc);

    void screen_output_end_info(const Simulation::DetQMC &dqmc);

    void screen_output_time();

    void screen_output_mpi(const std::string &master_proc_name, const int &world_size);

    void screen_output_params(const int &world_size, const Simulation::DetQMC &dqmc);

    void screen_output_observable(const Measure::MeasureData &obs, const std::string &obs_name);

    // void screen_output_observable_bin(const Measure::MeasureData &obs);

    // void screen_output_observable(const std::vector<Measure::MeasureData> &obs_vec);

    // void screen_output_observable_bin(const std::vector<Measure::MeasureData> &obs_vec);

} // namespace ScreenOutput

#endif // DQMC_HUBBARD_OUTPUT_H