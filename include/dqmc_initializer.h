#ifndef DQMC_INITIALIZER_H
#define DQMC_INITIALIZER_H
#pragma once


/**
  *  This header file defines the dqmc initializer class QuantumMonteCarlo::DqmcInitializer,
  *  which contains static member functions to initialize the dqmc modules altogether. 
  */

 #include <vector>

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
    using MomentumIndex = int;
    using MomentumIndexList = std::vector<int>;


    // ----------------------- Interface class QuantumMonteCarlo::DqmcInitializer ------------------------
    class DqmcInitializer {
        public:

            // set up measured momentum for momentum-dependent observables
            // note that the MomentumIndex should be provided by specific lattice module,
            // and do not mannully assign the input momentum to avoid unexpected mistakes.
            // some pre-designed interfaces of lattice momentum ( k stars ) are provided
            // in the derived Lattice classes for measuring usages.
            static void set_measured_momentum       ( MeasureHandler& meas_handler, 
                                                      const MomentumIndex& momentum_index );


            // set up list of measured momentum for momentum-dependent observables
            // the annotation above also applies to this function.
            static void set_measured_momentum_list  ( MeasureHandler& meas_handler, 
                                                      const MomentumIndexList& momentum_index_list );


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