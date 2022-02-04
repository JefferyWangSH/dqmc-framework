#include "output.h"
#include "detqmc.h"
#include "hubbard.h"

#include <iostream>
#include <fstream>

#include <boost/format.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time/local_time.hpp>


namespace FileOutput {

    void file_output_observable(const Measure::Observable<double> &obs, const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }
        
        // for specfic measuring object of double type, output mean value, error bar and relative error in order. 
        boost::format fmt_obs("%| 20.10f|%| 20.10f|%| 20.10f|");
        outfile << fmt_obs % obs.mean_value() % obs.error_bar() % (obs.error_bar()/obs.mean_value()) << std::endl;
        outfile.close();
    }

    void file_output_observable(const Measure::Observable<Eigen::VectorXd> &obs, const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }
        
        // output of vector type observable 
        boost::format fmt_info("%| 20d|");
        boost::format fmt_obs("%| 20d|%| 20.10f|%| 20.10f|%| 20.10f|");

        const int size = obs.mean_value().size();
        const Eigen::VectorXd relative_error = (obs.error_bar().array()/obs.mean_value().array()).matrix();
        outfile << fmt_info % size << std::endl;
        for (int i = 0; i < size; ++i) {
            outfile << fmt_obs % i % obs.mean_value()(i) % obs.error_bar()(i) % relative_error(i) << std::endl;
        }
        outfile.close();
    }

    void file_output_observable(const Measure::Observable<Eigen::MatrixXd> &obs, const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }
        
        // output of vector type observable 
        boost::format fmt_info("%| 20d|%| 20d|");
        boost::format fmt_obs("%| 20d|%| 20d|%| 20.10f|%| 20.10f|%| 20.10f|");

        const int row = obs.mean_value().rows();
        const int col = obs.mean_value().cols();
        const Eigen::MatrixXd relative_error = (obs.error_bar().array()/obs.mean_value().array()).matrix();
        outfile << fmt_info % row % col << std::endl;
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                outfile << fmt_obs % i % j % obs.mean_value()(i,j) % obs.error_bar()(i,j) % relative_error(i,j) << std::endl;
            }
        }
        outfile.close();
    }

    void file_output_observable_bin(const Measure::Observable<double> &obs, const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }

        // output bin data of observables
        boost::format fmt_info("%| 20d|");
        boost::format fmt_obs("%| 20d|%| 20.10f|");
        int nbin = obs.size_of_bin();
        outfile << fmt_info % nbin << std::endl;
        for (int bin = 0; bin < nbin; ++bin) {
            outfile << fmt_obs % bin % obs.bin_data()[bin] << std::endl;
        }
        outfile.close();
    }

    void file_output_observable_bin(const Measure::Observable<Eigen::VectorXd> &obs, const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }

        // output bin data of observables
        boost::format fmt_info("%| 20d|%| 20d|");
        boost::format fmt_obs("%| 20d|%| 20d|%| 20.10f|");
        int nbin = obs.size_of_bin();
        int size = obs.mean_value().size();
        outfile << fmt_info % size % nbin << std::endl;
        for (int i = 0; i < size; ++i) {
            for (int bin = 0; bin < nbin; ++bin) {
                outfile << fmt_obs % i % bin % obs.bin_data()[bin](i) << std::endl;
            }
        }
        outfile.close();
    }

    void file_output_observable_bin(const Measure::Observable<Eigen::MatrixXd> &obs, const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }

        // output bin data of observables
        boost::format fmt_info("%| 20d|%| 20d|%| 20d|");
        boost::format fmt_obs("%| 20d|%| 20d|%| 20d|%| 20.10f|");
        int nbin = obs.size_of_bin();
        int row = obs.mean_value().rows();
        int col = obs.mean_value().cols();
        outfile << fmt_info % row % col % nbin << std::endl;
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                for (int bin = 0; bin < nbin; ++bin) {
                    outfile << fmt_obs % i % j % bin % obs.bin_data()[bin](i,j) << std::endl;
                }
            }
        }
        outfile.close();
    }

    void file_output_tau(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }

        // output sequence of imaginary-time grids
        boost::format fmt_tau_info("%| 20d|%| 20.2f|");
        boost::format fmt_tau_seq("%| 20d|%| 20.10f|");
        outfile << fmt_tau_info % dqmc.hubbard->lt % dqmc.hubbard->beta << std::endl;
        for (int l = 0; l < dqmc.hubbard->lt; ++l) {
            outfile << fmt_tau_seq % l % (l*dqmc.hubbard->dtau) << std::endl;
        }
        outfile.close();
    }

    void file_output_aux_field(const Simulation::DetQMC &dqmc, const std::string &file_name, const int &mode) {
        std::ofstream outfile;
        if (mode) { outfile.open(file_name, std::ios::out | std::ios::app); }
        else { outfile.open(file_name, std::ios::out | std::ios::trunc); }

        // output current configuration of aux fields
        boost::format fmt_field_info("%| 20d|%| 20d|");
        boost::format fmt_field_seq("%| 20d|%| 20d|%| 20.1f|");
        outfile << fmt_field_info % dqmc.hubbard->lt % dqmc.hubbard->ls << std::endl;
        for (int l = 0; l < dqmc.hubbard->lt; ++l) {
            for (int i = 0; i < dqmc.hubbard->ls; ++i) {
                outfile << fmt_field_seq % l % i % (*dqmc.hubbard->s)(i, l) << std::endl;
            }
        }
        outfile.close();
    }

} // namespace FileOutput


namespace ScreenOutput {

    void screen_output_time() {
        // print current date and time
        auto current_time = boost::posix_time::second_clock::local_time();
        std::cout << " Current time : " << current_time << "\n" << std::endl;
    }

    void screen_output_mpi(const std::string &master_proc_name, const int &world_size) {
        // print information of processors
        std::cout << boost::format(" Distribute tasks to %s processors, with the master processor being %s. \n") % world_size % master_proc_name 
                  << std::endl;
    }

    void screen_output_params(const int &world_size, const Simulation::DetQMC &dqmc){
        // print simualtion parameters          
        boost::format fmt_param_int("%| 30s|%| 5s|%| 8d|");
        boost::format fmt_param_double("%| 30s|%| 5s|%| 8.3f|");
        boost::format fmt_param_k("%| 30s|%| 5s|%| 8.3f| pi, %.3f pi");
        const std::string joiner = "->";

        if (!dqmc.is_warm_up) { std::cout << " Configurations of aux fields read from input config file. \n" << std::endl;}
        else { std::cout << " Configurations of aux fields set to random. \n" << std::endl;}
        std::cout << " Initialization finished. \n\n"
                  << " The simulation is going to get started with parameters shown below : \n" << std::endl;
        std::cout << fmt_param_int % "Lattice length 'll'" % joiner % dqmc.hubbard->ll << std::endl;
        std::cout << fmt_param_int % "Imaginary-time length 'lt'" % joiner % dqmc.hubbard->lt << std::endl;
        std::cout << fmt_param_double % "Inverse temperature 'beta'" % joiner % dqmc.hubbard->beta << std::endl;
        std::cout << fmt_param_double % "Interaction strength 'U'" % joiner % dqmc.hubbard->u_int << std::endl;
        std::cout << fmt_param_double % "Chemical potential 'mu'" % joiner % dqmc.hubbard->mu << std::endl;
        std::cout << fmt_param_k % "Lattice momentum 'k'" % joiner % (dqmc.q[0]/M_PI) % (dqmc.q[1]/M_PI) << std::endl;
        std::cout << std::endl;

        std::cout << fmt_param_int % "Stablization pace 'nwrap'" % joiner % dqmc.nwrap << std::endl;
        std::cout << fmt_param_int % "Number of bins 'nbin'" % joiner % (dqmc.nbin * world_size) << std::endl;
        std::cout << fmt_param_int % "Sweeps per bin 'nsweep'" % joiner % dqmc.nsweep << std::endl << std::endl;
    }

    void screen_output_init_info(const std::string &master_proc_name, const int &world_size, const Simulation::DetQMC &dqmc) {
        // print current date and time
        screen_output_time();
        
        // print information of processors
        screen_output_mpi(master_proc_name, world_size);
        
        // print simualtion parameters 
        screen_output_params(world_size, dqmc);
    }

    void screen_output_end_info(const Simulation::DetQMC &dqmc) {
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(dqmc.end_t - dqmc.begin_t).count();
        const int day = std::floor((double)time / 86400000);
        const int hour = std::floor(((double)time/1000 - day * 86400) / 3600);
        const int minute = std::floor(((double)time/1000 - day * 86400 - hour * 3600) / 60);
        const double sec = (double)time/1000 - 86400 * day - 3600 * hour - 60 * minute;

        // print the time cost of simulation
        if ( day ) { std::cout << boost::format("\n The simulation finished in %d d %d h %d m %.2f s. \n") % day % hour % minute % sec << std::endl; }
        else if ( hour ) { std::cout << boost::format("\n The simulation finished in %d h %d m %.2f s. \n") % hour % minute % sec << std::endl; }
        else if ( minute ) { std::cout << boost::format("\n The simulation finished in %d m %.2f s. \n") % minute % sec << std::endl; }
        else { std::cout << boost::format("\n The simulation finished in %.2f s. \n") % sec << std::endl; }
        
        // print wrap error of the evaluation of Green's functions
        if (dqmc.is_eqtime_measure || dqmc.is_dynamic_measure || dqmc.is_warm_up) {
            std::cout << " Maximum of equal-time wrap error :  " << dqmc.hubbard->max_wrap_error_equal << std::endl;
            if (dqmc.is_dynamic_measure) {
                std::cout << "\n Maximum of dynamical wrap error :   " << dqmc.hubbard->max_wrap_error_dynamic << std::endl;
            }
            std::cout << std::endl;
        }
    }

    void screen_output_observable(const Measure::Observable<double> &obs, const std::string &obs_name) {
        boost::format fmt_obs("%| 30s|%| 5s|%| 17.12f| pm %.12f");
        const std::string joiner = "->";
        std::cout << fmt_obs % obs_name % joiner % obs.mean_value() % obs.error_bar() << std::endl;
    }

} // namespace ScreenOutput
