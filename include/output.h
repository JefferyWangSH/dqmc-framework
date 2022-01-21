#ifndef DQMC_HUBBARD_OUTPUT_H
#define DQMC_HUBBARD_OUTPUT_H
#pragma once

/**
  *  This head file includes declarations of subroutines for outputing DQMC measuring results,
  *  in form of outfile stream or terminal output.
  */

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/format.hpp>
#include "measure_data.h"


// forward declaration
namespace Simulation { class DetQMC; }

namespace FileOutput{

    template<typename DataStructure>
    void file_output_observable(const Measure::MeasureData<DataStructure> &obs, 
                                const std::string &file_name, const int &mode);

    template<typename DataStructure>
    void file_output_observable_bin(const Measure::MeasureData<DataStructure> &obs, 
                                    const std::string &file_name, const int &mode);

    template<typename DataStructure>
    void file_output_observable(const std::vector<Measure::MeasureData<DataStructure>> &obs_vec, 
                                const std::string &file_name, const int &mode);

    template<typename DataStructure>
    void file_output_observable_bin(const std::vector<Measure::MeasureData<DataStructure>> &obs_vec, 
                                    const std::string &file_name, const int &mode);

    void file_output_tau(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode);

    void file_output_aux_field(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode);



    /* Implements of template functions */

    template<typename DataStructure>
    void file_output_observable(const Measure::MeasureData<DataStructure> &obs, 
                                const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }
        
        // for a specfic measuring object, output mean value, error bar and relative error in order. 
        boost::format fmt_obs("%| 20.10f|%| 20.10f|%| 20.10f|");
        outfile << fmt_obs % obs.mean_value() % obs.error_bar() % (obs.error_bar()/obs.mean_value()) << std::endl;
        outfile.close();
    }

    template<typename DataStructure>
    void file_output_observable_bin(const Measure::MeasureData<DataStructure> &obs, 
                                    const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }

        // output bin data of observables
        boost::format fmt_bin_info("%| 20d|");
        boost::format fmt_bin_obs("%| 20d|%| 20.10f|");
        outfile << fmt_bin_info % obs.size_of_bin() << std::endl;
        for (int bin = 0; bin < obs.size_of_bin(); ++bin) {
            outfile << fmt_bin_obs % bin % obs.bin_data()[bin] << std::endl;
        }
        outfile.close();
    }

    template<typename DataStructure>
    void file_output_observable(const std::vector<Measure::MeasureData<DataStructure>> &obs_vec, 
                                const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }

        // for a list of observables, output mean value, error bar and relative error in order
        boost::format fmt_obs("%| 20d|%| 20.10f|%| 20.10f|%| 20.10f|");
        for (int i = 0; i < obs_vec.size(); ++i) {
            outfile << fmt_obs % i % obs_vec[i].mean_value() % obs_vec[i].error_bar() % (obs_vec[i].error_bar()/obs_vec[i].mean_value()) << std::endl;
        }
        outfile.close(); 
    }

    template<typename DataStructure>
    void file_output_observable_bin(const std::vector<Measure::MeasureData<DataStructure>> &obs_vec, 
                                    const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }

        // output bin data for a list of observables
        boost::format fmt_bin_info("%| 20d|%| 20d|");
        boost::format fmt_bin_obs("%| 20d|%| 20d|%| 20.10f|");
        outfile << fmt_bin_info % obs_vec.size() % obs_vec[0].size_of_bin() << std::endl;
        for (int i = 0; i < obs_vec.size(); ++i) {
            for (int bin = 0; bin < obs_vec[i].size_of_bin(); ++bin) {
                outfile << fmt_bin_obs % i % bin % obs_vec[i].bin_data()[bin] << std::endl;
            }
        }
        outfile.close();
    }

} // namespace FileOutput


namespace ScreenOutput {

    void screen_output_init_info(const std::string &master_proc_name, const int &world_size, const Simulation::DetQMC &dqmc);

    void screen_output_end_info(const Simulation::DetQMC &dqmc);

    void screen_output_time();

    void screen_output_mpi(const std::string &master_proc_name, const int &world_size);

    void screen_output_params(const int &world_size, const Simulation::DetQMC &dqmc);

    template<typename DataStructure>
    void screen_output_observable(const Measure::MeasureData<DataStructure> &obs, const std::string &obs_name);


    /* Implements of template functions */

    template<typename DataStructure>
    void screen_output_observable(const Measure::MeasureData<DataStructure> &obs, const std::string &obs_name) {
        boost::format fmt_obs("%| 30s|%| 5s|%| 17.12f| pm %.12f");
        const std::string joiner = "->";
        std::cout << fmt_obs % obs_name % joiner % obs.mean_value() % obs.error_bar() << std::endl;
    }

} // namespace ScreenOutput

#endif // DQMC_HUBBARD_OUTPUT_H