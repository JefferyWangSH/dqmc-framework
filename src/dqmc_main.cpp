#include "detqmc.h"
#include "random.h"
#include "output.h"
#include "measure_gather.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>

#include <mpi.h>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/program_options.hpp>
#include <boost/format.hpp>


/**
  *  TODO:
  *   1. get params from command lines, using boost (done)
  *   2. equal-time measurements of momentum distribution and spin-spin correlation (done)
  *   3. bin measurements (done)
  *   4. time-displaced green function and measurements (done)
  *   5. ******** Modify command console output ******** (done)
  *   6. attractive interaction U < 0 (done)
  *   7. determine the critical temperature of superconducting transition (done)
  *   8. reweighing for doped case (done)
  *   9. read aux field configurations from input file (done)
  *   10. checkerboard decomposition (done)
  *   11. openmp parallel sampling (missing)
  *   12. log output (done)
  *   13. run the program with bash scripts (done)
  *   14. MPI distributed programing (done)
  *   15. generalized model and lattice module supporting simulations of customized physical systems (missing)
  *   16. using fftw3 for fast fourier transformation of measurements in momentum space (missing)
  *   17. independent random (and seed) module (done)
  *   18. ...
  */


/**  The Main Program */
int main(int argc, char* argv[]) {

    /** model and controlling params */
    int ll = 4;
    int lt = 80;
    double beta = 4.0;
    double t = 1.0;
    double u = -4.0;
    double mu = 0.0;

    int nwrap = 10;
    int nwarm = (int)(4*ll*ll*beta);

    int nbin = 20;
    int nsweep = 100;
    int n_between_bins = 10;

    bool is_checkerboard = false;
    bool is_warm_up = true;
    bool is_eqtime_measure = true;
    bool is_dynamic_measure = true;

    std::string out_folder_path = "../results/example";

    std::vector<std::string> observable_supported = {"filling_number", "double_occupancy", "kinetic_energy",
                                                     "momentum_distribution", "local_spin_corr", 
                                                     "spin_density_structure_factor", "charge_density_structure_factor",
                                                     "s_wave_pairing_corr", 
                                                     "greens_functions", "density_of_states", 
                                                     "superfluid_stiffness", };
    std::vector<std::string> obs_list = observable_supported;

    /** read params from command line */
    boost::program_options::options_description opts("Program options");
    boost::program_options::variables_map vm;

    opts.add_options()
        ("help,h", "display this information")
        ("ll", boost::program_options::value<int>(&ll)->default_value(4), "spatial size of lattice, default: 4")
        ("lt", boost::program_options::value<int>(&lt)->default_value(80), "imaginary-time size of lattice, default: 80")
        ("beta", boost::program_options::value<double>(&beta)->default_value(4.0), "inverse temperature, default: 4.0")
        ("t", boost::program_options::value<double>(&t)->default_value(1.0), "hopping strength, default: 1.0")
        ("u", boost::program_options::value<double>(&u)->default_value(-4.0), "interaction strength, u > 0 for repulsive and u < 0 for attractive case, default: -4.0")
        ("mu", boost::program_options::value<double>(&mu)->default_value(0.0), "chemical potential, default: 0.0")
        ("checker-board", boost::program_options::value<bool>(&is_checkerboard)->default_value(false), "whether to perform checkerboard break-up, default: false")
        ("nwrap", boost::program_options::value<int>(&nwrap)->default_value(10), "pace of stabilization process, default: 10")
        ("nwarm", boost::program_options::value<int>(&nwarm)->default_value((int)(4*ll*ll*beta)), "number of warmup sweeps, default: 4*ll*ll*beta")
        ("nbin", boost::program_options::value<int>(&nbin)->default_value(20), "number of bins, default: 20")
        ("nsweep", boost::program_options::value<int>(&nsweep)->default_value(100), "number of measurement sweeps in a bin, default: 100")
        ("n-between-bins", boost::program_options::value<int>(&n_between_bins)->default_value(10), "number of sweeps between bins to avoid correlation, default: 10")
        ("warm-up", boost::program_options::value<bool>(&is_warm_up)->default_value(true), "whether to warm up, default: true")
        ("eqtime-measure", boost::program_options::value<bool>(&is_eqtime_measure)->default_value(true), "whether to do equal-time measurements, default: true")
        ("dynamic-measure", boost::program_options::value<bool>(&is_dynamic_measure)->default_value(true), "whether to do dynamic measurements, default: true")
        ("observable-list", boost::program_options::value<std::vector<std::string>>(&obs_list)->multitoken(), "list of physical observables to be measured, default: all")
        ("out-folder-path", boost::program_options::value<std::string>(&out_folder_path)->default_value("../results/example"), "path of the output folder, default: ../results/example");
    
    try {
        boost::program_options::store(parse_command_line(argc, argv, opts), vm);
    }
    catch (...) {
        std::cerr << " Undefined options got from command line!\n"<< std::endl;
        exit(1);
    }
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        std::cerr << argv[0] << std::endl;
        std::cerr << opts << std::endl;
        return 0;
    }
    if (vm.count("observable-list")) {
        if (obs_list.size() == 1 && obs_list[0] == "all") {
            obs_list = observable_supported;
        }
        if (obs_list.size() == 1 && obs_list[0] == "none") {
            obs_list = std::vector<std::string>();
        }
    }
    if ((!vm["ll"].defaulted() || !vm["beta"].defaulted()) && vm["nwarm"].defaulted()) {
        nwarm = 4 * vm["ll"].as<int>() * vm["ll"].as<int>() * (int)vm["beta"].as<double>();
    }


    /** DQMC Simulation */
    Simulation::DetQMC *dqmc;
    dqmc = new Simulation::DetQMC();

    /** distributed parallelization using MPI */
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    const int master = 0;
    const int rank = world.rank();

    // arrange bins measurements over processors
    const int bins_per_proc = (nbin % world.size() == 0)? nbin/world.size() : nbin/world.size()+1;

    // set up unique random seed for each processor
    Random::set_seed(rank);

    // usage example
    dqmc->set_model_params(ll, lt, beta, t, u, mu, nwrap);
    dqmc->set_Monte_Carlo_params(nwarm, bins_per_proc, nsweep, n_between_bins);
    dqmc->set_controlling_params(is_warm_up, is_eqtime_measure, is_dynamic_measure, is_checkerboard);
    dqmc->set_observable_list(obs_list);
    dqmc->set_lattice_momentum((Eigen::VectorXd(2) << 1.0*M_PI, 1.0*M_PI).finished());

    // initialize output folder
    if (rank == master) {
        if ( access(out_folder_path.c_str(), 0) != 0 ) {
            const std::string command = "mkdir " + out_folder_path;
            if ( system(command.c_str()) != 0 ) {
                std::cerr << boost::format(" Fail to creat folder at %s . \n") % out_folder_path << std::endl;
            }
        }
    }
    // read in configurations of aux fields
    if (!is_warm_up) {    
        dqmc->set_aux_field_configs(out_folder_path + "/config.dat");
    }

    // initialization
    dqmc->initial();
    
    // the master processor
    if (rank == master) {
        // output the information of simulation
        ScreenOutput::screen_output_init_info(env.processor_name(), world.size(), *dqmc);
    }

    // MC simulation process ...
    const bool show_running_process = (rank == master);
    dqmc->run(show_running_process);

    // output the ending information, include time cost and wrap errors
    if (rank == master) {
        ScreenOutput::screen_output_end_info(*dqmc);
    }

    // analyse statistics
    dqmc->analyse_stats();

    // collect data and output results
    if (dqmc->measure) {
        // collect data from all processors
        // Todo: redesign interface of measure and container class
        Measure::GatherMPI gather_mpi{};
        gather_mpi.gather_from_all_processors(world, *dqmc->measure);

        if (rank == master) {
            // display measuring results on terminal
            // print results of all double-type observables
            // equal-time measurements
            if (dqmc->measure->find("filling_number")) {
                ScreenOutput::screen_output_observable(dqmc->measure->find_double_obs("filling_number"), "Filling number");
            }
            if (dqmc->measure->find("double_occupancy")) {
                ScreenOutput::screen_output_observable(dqmc->measure->find_double_obs("double_occupancy"), "Double occupancy");
            }
            if (dqmc->measure->find("kinetic_energy")) {
                ScreenOutput::screen_output_observable(dqmc->measure->find_double_obs("kinetic_energy"), "Kinetic energy");
            }
            if (dqmc->measure->find("momentum_distribution")) {
                ScreenOutput::screen_output_observable(dqmc->measure->find_double_obs("momentum_distribution"), "Momentum distribution");
            }
            if (dqmc->measure->find("local_spin_corr")) {
                ScreenOutput::screen_output_observable(dqmc->measure->find_double_obs("local_spin_corr"), "Local spin correlation");
            }
            if (dqmc->measure->find("spin_density_structure_factor")) {
                ScreenOutput::screen_output_observable(dqmc->measure->find_double_obs("spin_density_structure_factor"), "SDW order parameter");
            }
            if (dqmc->measure->find("charge_density_structure_factor")) {
                ScreenOutput::screen_output_observable(dqmc->measure->find_double_obs("charge_density_structure_factor"), "CDW order parameter");
            }           
            if (dqmc->measure->find("s_wave_pairing_corr")) {
                ScreenOutput::screen_output_observable(dqmc->measure->find_double_obs("s_wave_pairing_corr"), "S-wave pairing correlation");
            } 
            if (dqmc->measure->is_eqtime_measure()) {
                ScreenOutput::screen_output_observable(dqmc->measure->sign_eqtime(), "Averaged sign (equal-time)");
            }

            // dynamical measurements
            if (dqmc->measure->find("superfluid_stiffness")) {
                ScreenOutput::screen_output_observable(dqmc->measure->find_double_obs("superfluid_stiffness"), "Superfluid stiffness");
            }
            if (dqmc->measure->is_dynamic_measure()) {
                ScreenOutput::screen_output_observable(dqmc->measure->sign_dynamic(), "Averaged sign (dynamic)");
            }

            // file output of measuring results
            // customize here
            if (dqmc->measure->find("s_wave_pairing_corr")) {
                FileOutput::file_output_observable(dqmc->measure->find_double_obs("s_wave_pairing_corr"), out_folder_path + "/swave.dat", 0);
                FileOutput::file_output_observable_bin(dqmc->measure->find_double_obs("s_wave_pairing_corr"), out_folder_path + "/swave.bin", 0);
            }
            if (dqmc->measure->find("superfluid_stiffness")) {
                FileOutput::file_output_observable(dqmc->measure->find_double_obs("superfluid_stiffness"), out_folder_path + "/sf.dat", 0);
                FileOutput::file_output_observable_bin(dqmc->measure->find_double_obs("superfluid_stiffness"), out_folder_path + "/sf.bin", 0);
            }
            if (dqmc->measure->find("greens_functions")) {
                FileOutput::file_output_observable(dqmc->measure->find_vector_obs("greens_functions"), out_folder_path + "/greens.dat", 0);
                FileOutput::file_output_observable_bin(dqmc->measure->find_vector_obs("greens_functions"), out_folder_path + "/greens.bin", 0);
            }

        }
    }


    // configuration output of aux fields
    if (rank == master) {
        FileOutput::file_output_aux_field(*dqmc, out_folder_path + "/config.dat", 0);
    }

    // // file output of imaginary-time grids
    // if (rank == master) {
    //     FileOutput::file_output_tau(*dqmc, out_folder_path + "/tau.dat", 0)
    // }

    delete dqmc;
    return 0;
}
