#include "dqmc_initializer.h"
#include "dqmc_walker.h"
#include "svd_stack.h"

#include "model/model_base.h"
#include "model/repulsive_hubbard.h"
#include "model/attractive_hubbard.h"

#include "lattice/lattice_base.h"
#include "lattice/square.h"
#include "lattice/cubic.h"
#include "lattice/honeycomb.h"

#include "checkerboard/checkerboard_base.h"
#include "checkerboard/square.h"
#include "checkerboard/cubic.h"

#include "measure/measure_handler.h"

#include "utils/toml.hpp"
#include <iostream>


namespace QuantumMonteCarlo {


    void DqmcInitializer::parse_toml_config(  std::string_view toml_config,
                                              int world_size,
                                              ModelBasePtr& model, 
                                              LatticeBasePtr& lattice, 
                                              DqmcWalkerPtr& walker,
                                              MeasureHandlerPtr& meas_handler,
                                              CheckerBoardBasePtr& checkerboard )
    {
        // parse the configuration file
        auto config = toml::parse_file( toml_config );


        // --------------------------------------------------------------------------------------------------
        //                                      Parse the Model module
        // --------------------------------------------------------------------------------------------------
        // create the model object and set up parameters case by case
        const std::string_view model_type = config["Model"]["type"].value_or("RepulsiveHubbard");
        
        // -----------------------------------  Repulsive Hubbard model  ------------------------------------
        if ( model_type == "RepulsiveHubbard" ) 
        {
            const double hopping_t          = config["Model"]["Params"]["hopping_t"].value_or(1.0);
            const double onsite_u           = config["Model"]["Params"]["onsite_u"].value_or(4.0);
            const double chemical_potential = config["Model"]["Params"]["chemical_potential"].value_or(0.0);
            
            if ( model ) { model.reset(); }
            model = std::make_unique<Model::RepulsiveHubbard>();
            model->set_model_params( hopping_t, onsite_u, chemical_potential );

        }

        // -----------------------------------  Attractive Hubbard model  -----------------------------------
        else if ( model_type == "AttractiveHubbard" ) 
        {
            const double hopping_t          = config["Model"]["Params"]["hopping_t"].value_or(1.0);
            const double onsite_u           = config["Model"]["Params"]["onsite_u"].value_or(4.0);
            const double chemical_potential = config["Model"]["Params"]["chemical_potential"].value_or(0.0);
            
            if ( model ) { model.reset(); }
            model = std::make_unique<Model::AttractiveHubbard>();
            model->set_model_params( hopping_t, onsite_u, chemical_potential );
        }

        else
        {
            std::cerr << "QuantumMonteCarlo::DqmcInitializer::parse_toml_config(): "
                      << "undefined model type \'" << model_type << "\', please check the config." << std::endl; 
            exit(1);
        }


        // --------------------------------------------------------------------------------------------------
        //                                    Parse the Lattice module
        // --------------------------------------------------------------------------------------------------
        // create the lattice object and set up parameters
        const std::string_view lattice_type = config["Lattice"]["type"].value_or("Square");

        // -------------------------------------  2D Square lattice  ----------------------------------------
        if ( lattice_type == "Square" ) 
        {   
            // parse size of the lattice
            std::vector<int> lattice_size;
            toml::array* lattice_arr = config["Lattice"]["cell"].as_array();
            if ( lattice_arr && lattice_arr->size() == 2 && lattice_arr->is_homogeneous<int64_t>() ) {
                lattice_size.reserve(lattice_arr->size());
                for ( auto&& el : *lattice_arr ) {
                    if ( el.value_or(0) < 1 ) {
                        std::cerr << "QuantumMonteCarlo::DqmcInitializer::parse_toml_config(): "
                                  << "the input cell of 2d square lattice should be a vector containing two positive intergers, "
                                  << "please check the config." << std::endl;
                        exit(1);
                    }
                    lattice_size.emplace_back(el.value_or(0));
                }
            }
            else {
                std::cerr << "QuantumMonteCarlo::DqmcInitializer::parse_toml_config(): "
                          << "the input cell of 2d square lattice should be a vector containing two positive intergers, "
                          << "please check the config." << std::endl;
                exit(1);
            }

            // create 2d square lattice object
            if ( lattice ) { lattice.reset(); }
            lattice = std::make_unique<Lattice::Square>();
            lattice->set_lattice_params( lattice_size );

            // initial lattice module in place
            if ( !lattice->InitialStatus() ) { lattice->initial(); }
        }

        // -------------------------------------  3D Cubic lattice  -----------------------------------------
        else if ( lattice_type == "Cubic" )
        {
            // todo
        }

        // ------------------------------------  2D Honeycomb lattice  --------------------------------------
        else if ( lattice_type == "Honeycomb" )
        {
            // todo
        }

        else 
        {
            std::cerr << "QuantumMonteCarlo::DqmcInitializer::parse_toml_config(): "
                      << "undefined lattice type \'" << lattice_type << "\', please check the config." << std::endl; 
            exit(1);
        }


        // --------------------------------------------------------------------------------------------------
        //                                  Parse the CheckerBoard module
        // --------------------------------------------------------------------------------------------------
        // note that the checkerboard method is currently only implemented for 2d square lattice
        const bool is_checker_board = config["CheckerBoard"]["whether_or_not"].value_or(false);
        
        if ( checkerboard ) { checkerboard.reset(); }
        if ( is_checker_board ) {
            if ( lattice_type == "Square" ) {
                checkerboard = std::make_unique<CheckerBoard::Square>();
            }
            else {
                std::cerr << "QuantumMonteCarlo::DqmcInitializer::parse_toml_config(): "
                          << "the checkerboard method is currently only implemented for 2d square lattice, "
                          << "please check the config." << std::endl;
                exit(1);
            }
        }


        // --------------------------------------------------------------------------------------------------
        //                                   Parse the DqmcWalker module
        // --------------------------------------------------------------------------------------------------
        const int beta = config["MonteCarlo"]["beta"].value_or(4.0);
        const double time_size = config["MonteCarlo"]["time_size"].value_or(80);
        const int stabilization_pace = config["MonteCarlo"]["stabilization_pace"].value_or(10);

        // create dqmc walker and set up parameters
        if ( walker ) { walker.reset(); }
        walker = std::make_unique<DqmcWalker>();
        walker->set_physical_params( beta, time_size );
        walker->set_stabilization_pace( stabilization_pace );


        // --------------------------------------------------------------------------------------------------
        //                                Parse the Measure Handler module
        // --------------------------------------------------------------------------------------------------
        const int sweeps_warmup = config["Measure"]["sweeps_warmup"].value_or(512);
        const int bin_num = config["Measure"]["bin_num"].value_or(20);
        const int bin_size = config["Measure"]["bin_size"].value_or(100);
        const int sweeps_between_bins = config["Measure"]["sweeps_between_bins"].value_or(20);
        
        // parse obervable lists
        std::vector<std::string> observables;
        toml::array* observable_arr = config["Measure"]["observables"].as_array();
        if ( observable_arr && observable_arr->is_homogeneous<std::string>() ) {
            observables.reserve(observable_arr->size());
            for ( auto&& el : *observable_arr ) {
                observables.emplace_back(el.value_or(""));
            }
        }
        else {
            std::cerr << "QuantumMonteCarlo::DqmcInitializer::parse_toml_config(): "
                      << "undefined observables, please check the config." << std::endl;
            exit(1);
        }

        // deal with special keywords ( all/All , none/None )
        if ( observables.size() == 1 ) {
            if ( observables[0] == "all" || observables[0] == "All" ) {
                observables = Measure::MeasureHandler::ObservableAll;
            }
            else if ( observables[0] == "none" || observables[0] == "None" ) {
                observables = {};
            }
        }

        // special observables, e.g. superfluid stiffness, are only supported for specific lattice type.
        if ( lattice_type != "Square" ) { 
            observables.erase( std::remove( std::begin(observables), std::end(observables), "superfluid_stiffness" ), 
                               std::end(observables) ); 
        }
        
        // create measure handler and set up parameters
        if ( meas_handler ) { meas_handler.reset(); }
        meas_handler = std::make_unique<Measure::MeasureHandler>();

        // send measuring tasks to a set of processes
        const int bins_per_proc = (bin_num % world_size == 0)? bin_num/world_size : bin_num/world_size+1;
        meas_handler->set_measure_params( sweeps_warmup, bins_per_proc, bin_size, sweeps_between_bins );
        meas_handler->set_observables( observables );


        // --------------------------------------------------------------------------------------------------
        //                                Parse the input Momentum parmas
        // --------------------------------------------------------------------------------------------------
        // set up momentum and momentum list for measurements
        const std::string_view momentum = config["Lattice"]["momentum"].value_or("");
        const std::string_view momentum_list = config["Lattice"]["momentum_list"].value_or("");
        
        // make sure that the lattice module is initialized ahead
        if ( lattice->InitialStatus() ) 
        {   
            // -----------------------------------  2D Square lattice  --------------------------------------
            if ( lattice_type == "Square" ) {
                // covert base class pointer to that of the derived square lattice class
                if ( const auto square_lattice = dynamic_cast<const Lattice::Square*>(lattice.get()) ) {

                    if ( momentum == "GammaPoint" )  { meas_handler->set_measured_momentum( square_lattice->GammaPointIndex() ); }
                    else if ( momentum == "MPoint" ) { meas_handler->set_measured_momentum( square_lattice->MPointIndex() ); }
                    else if ( momentum == "XPoint" ) { meas_handler->set_measured_momentum( square_lattice->XPointIndex() ); }
                    else { 
                        std::cerr << "QuantumMonteCarlo::DqmcInitializer::parse_toml_config(): "
                                  << "undefined momentum \'" << momentum << "\' for 2d square lattice, "
                                  << "please check the config." << std::endl;
                        exit(1);
                    }

                    if ( momentum_list == "KstarsAll" ) { meas_handler->set_measured_momentum_list( square_lattice->kStarsIndex() ); }
                    else if ( momentum_list == "DeltaLine" ) { meas_handler->set_measured_momentum_list( square_lattice->DeltaLineIndex() ); }
                    else if ( momentum_list == "ZLine" ) { meas_handler->set_measured_momentum_list( square_lattice->ZLineIndex() ); }
                    else if ( momentum_list == "SigmaLine" ) { meas_handler->set_measured_momentum_list( square_lattice->SigmaLineIndex() ); }
                    else if ( momentum_list == "Gamma2X2M2GammaLoop" ) { meas_handler->set_measured_momentum_list( square_lattice->Gamma2X2M2GammaLoopIndex() ); }
                    else { 
                        std::cerr << "QuantumMonteCarlo::DqmcInitializer::parse_toml_config(): "
                                  << "undefined momentum list \'" << momentum_list << "\' for 2d square lattice, "
                                  << "please check the config." << std::endl;
                        exit(1);
                    }
                }
                else {
                    std::cerr << "QuantumMonteCarlo::DqmcInitializer::parse_toml_config(): "
                              << "fail to convert \'Lattice::LatticeBase\' to \'Lattice::Square\'." << std::endl;
                    exit(1);
                }

            }

            // -----------------------------------  3D Cubic lattice  ---------------------------------------
            if ( lattice_type == "Cubic" ) {
                // todo
            }

            // ----------------------------------  2D Honeycomb lattice  ------------------------------------
            if ( lattice_type == "Honeycomb" ) {
                // todo
            }

        }
        

    }


    void DqmcInitializer::initial_modules( ModelBase& model, 
                                           LatticeBase& lattice, 
                                           DqmcWalker& walker,
                                           MeasureHandler& meas_handler )
    {
        // make sure that the module objects have been created,
        // and the parameters are setup correctly in advance.
        // notice that the orders of initializations below are important.

        // initialize lattice module
        if ( !lattice.InitialStatus() ) { lattice.initial(); }

        // initialize MeasureHandler module
        meas_handler.initial( lattice, walker );

        // initialize dqmcWalker module
        walker.initial( lattice, meas_handler );

        // initialize model module
        // naively link
        model.initial( lattice, walker );
        model.link();
    }


    void DqmcInitializer::initial_modules( ModelBase& model, 
                                           LatticeBase& lattice, 
                                           DqmcWalker& walker,
                                           MeasureHandler& meas_handler,
                                           CheckerBoardBase& checkerboard )
    {
        // make sure that the module objects have been created,
        // and the parameters are setup correctly in advance.
        // notice that the orders of initializations below are important.

        // initialize lattice module
        if ( !lattice.InitialStatus() ) { lattice.initial(); }

        // initialize MeasureHandler module
        meas_handler.initial( lattice, walker );

        // initialize dqmcWalker module
        walker.initial( lattice, meas_handler );

        // initialize model module
        model.initial( lattice, walker );

        // initialize checkerboard module and link to the model class
        checkerboard.set_checkerboard_params( lattice, model, walker );
        checkerboard.initial();
        model.link( checkerboard );
    }


    void DqmcInitializer::initial_dqmc( ModelBase& model, 
                                        LatticeBase& lattice, 
                                        DqmcWalker& walker,
                                        MeasureHandler& meas_handler )
    {
        // this subroutine should be called after the initial 
        // configuration of the bosonic fields have been determined, 
        // either randomly initialized or read from a input config file.
        // SvdStack class are initialized and the greens functions 
        // for the initial bosonic fields are computed in this function.
        walker.initial_svd_stacks( lattice, model );
        walker.initial_greens_functions();
        walker.initial_config_sign();
    }


} // namespace QuantumMonteCarlo