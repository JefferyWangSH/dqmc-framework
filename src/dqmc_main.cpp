#include <memory>
#include <string>
#include <iostream>
#include <fstream>

#include <mpi.h>
#include <boost/mpi.hpp>
#include <boost/format.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/program_options.hpp>

#include "model/model_base.h"
#include "lattice/lattice_base.h"
#include "checkerboard/checkerboard_base.h"
#include "measure/measure_handler.h"
#include "dqmc.h"
#include "dqmc_walker.h"

#include "dqmc_initializer.h"
#include "dqmc_io.h"
#include "svd_stack.h"
#include "random.h"
#include "utils/mpi.hpp"



// the main program
int main( int argc, char* argv[] ) {
    
    // ------------------------------------------------------------------------------------------------
    //                                 Initialize MPI environment
    // ------------------------------------------------------------------------------------------------
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;
    const int master = 0;
    const int rank = world.rank();


    // ------------------------------------------------------------------------------------------------
    //                                      Program options
    // ------------------------------------------------------------------------------------------------
    
    std::string config_file;
    std::string fields_file;
    std::string out_path;
    
    // read parameters from the command line 
    boost::program_options::options_description opts("Program options");
    boost::program_options::variables_map vm;

    opts.add_options()
        (   "help,h", "display this information" )
        (   "config,c", 
            boost::program_options::value<std::string>(&config_file)->default_value("../example/config.toml"), 
            "path of the configuration file, default: ../example/config.toml" )
        (   "output,o",
            boost::program_options::value<std::string>(&out_path)->default_value("../example"), 
            "folder path which stores the output of measuring results, default: ../example" )
        (   "fields,f",
            boost::program_options::value<std::string>(&fields_file), 
            "path of the configurations of auxiliary fields, if not assigned the fields are to be set randomly." );
    
    // parse the command line options
    try {
        boost::program_options::store(parse_command_line(argc, argv, opts), vm);
    }
    catch (...) {
        std::cerr << " main(): undefined options got from command line." << std::endl; exit(1);
    }
    boost::program_options::notify(vm);

    // show the helping messages
    if (vm.count("help")) {
        std::cerr << argv[0] << "\n" << opts << std::endl;
        return 0;
    }

    // initialize the output folder, create if not exist
    if ( rank == master ) {
        if ( access(out_path.c_str(), 0) != 0 ) {
            const std::string command = "mkdir -p " + out_path;
            if ( system(command.c_str()) != 0 ) {
                std::cerr << boost::format(" main(): fail to creat folder at %s .\n") % out_path 
                          << std::endl;
                exit(1);
            }
        }
    }


    // ------------------------------------------------------------------------------------------------
    //                                Output current date and time
    // ------------------------------------------------------------------------------------------------
    if ( rank == master ) {
        const auto current_time = boost::posix_time::second_clock::local_time();
        std::cout << boost::format(" Current time: %s \n") % current_time << std::endl;
    }


    // ------------------------------------------------------------------------------------------------
    //                                Output MPI and hardware info
    // ------------------------------------------------------------------------------------------------
    if ( rank == master ) {
        // print MPI and hardware information
        boost::format fmt_mpi(" Distribute tasks to %s processes, with the master process being %s. \n");
        std::cout << fmt_mpi % world.size() % env.processor_name() << std::endl;
    }


    // ------------------------------------------------------------------------------------------------
    //                                 Process of DQMC simulation
    // ------------------------------------------------------------------------------------------------

    // set up random seeds for different processes
    Utils::Random::set_seed( std::time(nullptr) + rank );
    // // fixed random seed for debug
    // Utils::Random::set_seed( 12345 );
    
    
    // -----------------------------------  Initializations  ------------------------------------------

    // create dqmc module objects
    std::unique_ptr<Model::ModelBase> model;
    std::unique_ptr<Lattice::LatticeBase> lattice;
    std::unique_ptr<QuantumMonteCarlo::DqmcWalker> walker;
    std::unique_ptr<Measure::MeasureHandler> meas_handler;
    std::unique_ptr<CheckerBoard::CheckerBoardBase> checkerboard;

    // parse parmas from the configuation file
    QuantumMonteCarlo::DqmcInitializer::parse_toml_config
        ( 
            config_file, world.size(),
            model, lattice, walker, meas_handler, checkerboard 
        );

    // initialize modules
    if ( checkerboard ) { 
        // using checkerboard break-up
        QuantumMonteCarlo::DqmcInitializer::initial_modules( *model, *lattice, *walker, *meas_handler, *checkerboard ); 
    }
    else { 
        // without checkerboard break-up
        QuantumMonteCarlo::DqmcInitializer::initial_modules( *model, *lattice, *walker, *meas_handler ); 
    }

    if ( fields_file.empty() ) {
        // randomly initialize the auxiliary fields with no input fields configs
        model->set_bosonic_fields_to_random();
        if ( rank == master ) { 
            std::cout << " Configurations of auxiliary fields set to random. \n" << std::endl; 
        }
    }
    else {
        QuantumMonteCarlo::DqmcIO::read_bosonic_fields_from_file( fields_file, *model);
        if ( rank == master ) { 
            std::cout << " Configurations of auxiliary fields read from input config file. \n" << std::endl; 
        }
    }

    // initialize dqmc, preparing for the simulation
    QuantumMonteCarlo::DqmcInitializer::initial_dqmc( *model, *lattice, *walker, *meas_handler );

    if ( rank == master ) {
        std::cout << " Initialization finished. \n\n" 
                  << " The simulation is going to get started with parameters shown below : \n"
                  << std::endl;
    }

    // output the initialization info
    if ( rank == master ) {
        QuantumMonteCarlo::DqmcIO::output_init_info 
            ( 
                std::cout, world.size(), 
                *model, *lattice, *walker, *meas_handler, checkerboard 
            );
    }

    // set up progress bar
    QuantumMonteCarlo::Dqmc::show_progress_bar( (rank == master) );
    QuantumMonteCarlo::Dqmc::progress_bar_format( 60, '=', ' ' );
    QuantumMonteCarlo::Dqmc::set_refresh_rate( 10 );


    // ---------------------------------  Crucial simulation steps  ------------------------------------

    // the dqmc simulation start
    QuantumMonteCarlo::Dqmc::timer_begin();
    QuantumMonteCarlo::Dqmc::thermalize( *walker, *model, *lattice, *meas_handler );
    QuantumMonteCarlo::Dqmc::measure( *walker, *model, *lattice, *meas_handler );

    // gather observable objects from other processes
    Utils::MPI::mpi_gather( world, *meas_handler );

    // perform the analysis
    QuantumMonteCarlo::Dqmc::analyse( *meas_handler );

    // end the timer
    QuantumMonteCarlo::Dqmc::timer_end();

    // output the ending info
    if ( rank == master ) {
        QuantumMonteCarlo::DqmcIO::output_ending_info( std::cout, *walker );
    }


    // ---------------------------------  Output measuring results  ------------------------------------

    // screen output the results of scalar observables
    if ( rank == master ) 
    {
        if ( meas_handler->find("equaltime_sign") ) {
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("equaltime_sign") );
        }

        if ( meas_handler->find("dynamic_sign") ) {
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("dynamic_sign") );
        }

        std::cout << std::endl;
      
        if ( meas_handler->find("filling_number") ) {
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("filling_number") );
        }

        if ( meas_handler->find("double_occupancy") ) {
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("double_occupancy") );
        }

        if ( meas_handler->find("kinetic_energy") ) {
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("kinetic_energy") );
        }

        if ( meas_handler->find("local_spin_corr") ) {
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("local_spin_corr") );
        }

        if ( meas_handler->find("momentum_distribution") ) {
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("momentum_distribution") );
        }

        if ( meas_handler->find("spin_density_structure_factor") ) {
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("spin_density_structure_factor") );
        }

        if ( meas_handler->find("charge_density_structure_factor") ) {
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("charge_density_structure_factor") );
        }

        if ( meas_handler->find("s_wave_pairing_corr") ) {
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("s_wave_pairing_corr") );
        }

        if ( meas_handler->find("superfluid_stiffness") ) {  
            QuantumMonteCarlo::DqmcIO::output_observable( 
                std::cout, meas_handler->find<Observable::ScalarObs>("superfluid_stiffness") );
        }
    }


    // file output 
    if ( rank == master ) {
        
        std::ofstream outfile;

        // output the auxiliary fields configurations
        // if there exist input file of fields configs, overwrite it.
        // otherwise the field configs are stored under the output folder.
        const auto fields_out = ( fields_file.empty() )? out_path + "/fields.out" : fields_file;
        outfile.open(fields_out, std::ios::trunc);
        QuantumMonteCarlo::DqmcIO::output_bosonic_fields( outfile, *model );
        outfile.close();

        // output the k stars
        outfile.open(out_path + "/kstars.out", std::ios::trunc);
        QuantumMonteCarlo::DqmcIO::output_k_stars( outfile, *lattice );
        outfile.close();

        // output the imaginary-time grids
        outfile.open(out_path + "/tgrids.out", std::ios::trunc);
        QuantumMonteCarlo::DqmcIO::output_imaginary_time_grids( outfile, *walker );
        outfile.close();

        // output measuring results of the observables

        // s wave pairing correlation functions
        if ( meas_handler->find("s_wave_pairing_corr") ) {
            // output of means and errors
            outfile.open(out_path + "/swave.out", std::ios::trunc);
            QuantumMonteCarlo::DqmcIO::output_observable(
                outfile, meas_handler->find<Observable::ScalarObs>("s_wave_pairing_corr") );
            outfile.close();

            // output of raw data in terms of bins
            outfile.open(out_path + "/swave.bins.out", std::ios::trunc);
            QuantumMonteCarlo::DqmcIO::output_observable_in_bins(
                outfile, meas_handler->find<Observable::ScalarObs>("s_wave_pairing_corr") );
            outfile.close();
        }

        // density of states
        if ( meas_handler->find("density_of_states") ) {
            // output of means and errors
            outfile.open(out_path + "/dos.out", std::ios::trunc);
            QuantumMonteCarlo::DqmcIO::output_observable(
                outfile, meas_handler->find<Observable::VectorObs>("density_of_states") );
            outfile.close();

            // output of raw data in terms of bins
            outfile.open(out_path + "/dos.bins.out", std::ios::trunc);
            QuantumMonteCarlo::DqmcIO::output_observable_in_bins(
                outfile, meas_handler->find<Observable::VectorObs>("density_of_states") );
            outfile.close();
        }

        // dynamical green's function in the reciprocal space
        if ( meas_handler->find("greens_functions") ) {
            // output of means and errors
            outfile.open(out_path + "/greens.out", std::ios::trunc);
            QuantumMonteCarlo::DqmcIO::output_observable(
                outfile, meas_handler->find<Observable::MatrixObs>("greens_functions") );
            outfile.close();

            // output of raw data in terms of bins
            outfile.open(out_path + "/greens.bins.out", std::ios::trunc);
            QuantumMonteCarlo::DqmcIO::output_observable_in_bins(
                outfile, meas_handler->find<Observable::MatrixObs>("greens_functions") );
            outfile.close();
        }

    }



    return 0;

}
