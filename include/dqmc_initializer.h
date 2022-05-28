#ifndef DQMC_INITIALIZER_H
#define DQMC_INITIALIZER_H
#pragma once


/**
  *  This header file defines the dqmc initializer class QuantumMonteCarlo::DqmcInitializer,
  *  which contains static member functions to initialize the dqmc modules altogether. 
  */

 #include <vector>
 #include <string_view>
 #include <memory>

namespace Lattice { class LatticeBase; }
namespace Model { class ModelBase; }
namespace Measure { class MeasureHandler; }
namespace CheckerBoard { class CheckerBoardBase; }


namespace QuantumMonteCarlo {

    // forward declaration
    class DqmcWalker;

    using LatticeBase = Lattice::LatticeBase;
    using ModelBase = Model::ModelBase;
    using CheckerBoardBase = CheckerBoard::CheckerBoardBase;
    using MeasureHandler = Measure::MeasureHandler;
    
    using LatticeBasePtr = std::unique_ptr<Lattice::LatticeBase>;
    using ModelBasePtr = std::unique_ptr<Model::ModelBase>;
    using CheckerBoardBasePtr = std::unique_ptr<CheckerBoard::CheckerBoardBase>;
    using MeasureHandlerPtr = std::unique_ptr<Measure::MeasureHandler>;
    using DqmcWalkerPtr = std::unique_ptr<DqmcWalker>;
    
    using MomentumIndex = int;
    using MomentumIndexList = std::vector<int>;


    // ----------------------- Interface class QuantumMonteCarlo::DqmcInitializer ------------------------
    class DqmcInitializer {
        public:

            // parse parmameters from the toml configuration file
            // create modules and setup module parameters according to the input configurations
            static void parse_toml_config           ( std::string_view toml_config,
                                                      LatticeBasePtr& lattice, 
                                                      ModelBasePtr& model, 
                                                      DqmcWalkerPtr& walker,
                                                      MeasureHandlerPtr& meas_handler,
                                                      CheckerBoardBasePtr& checkerboard );


            // initialize modules including Lattice, Model, DqmcWalker and MeasureHandler
            // without checkerboard breakups.
            static void initial_modules             ( LatticeBase& lattice, 
                                                      ModelBase& model, 
                                                      DqmcWalker& walker,
                                                      MeasureHandler& meas_handler );
            

            // initialize modules including Lattice, Model, DqmcWalker and MeasureHandler
            // with checkerboard breakups. 
            static void initial_modules             ( LatticeBase& lattice, 
                                                      ModelBase& model, 
                                                      DqmcWalker& walker,
                                                      MeasureHandler& meas_handler,
                                                      CheckerBoardBase& checkerboard );


            // prepare for the dqmc simulation,
            // especially initializing the greens functions and SVD stacks
            static void initial_dqmc                ( LatticeBase& lattice, 
                                                      ModelBase& model,
                                                      DqmcWalker& walker,
                                                      MeasureHandler& meas_handler );

    };

} // namespace QuantumMonteCarlo


#endif // DQMC_INITIALIZER_H