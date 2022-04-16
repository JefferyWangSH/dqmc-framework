#ifndef DQMC_INITIALIZER_H
#define DQMC_INITIALIZER_H
#pragma once


/**
  *  This header file defines the dqmc initializer class QuantumMonteCarlo::DqmcInitializer,
  *  which contains static member functions to initialize the dqmc modules altogether. 
  */


namespace Lattice { class LatticeBase; }
namespace Model { class ModelBase; }
namespace Measure { class MeasureHandler; }


namespace QuantumMonteCarlo {

    // forward declaration
    class DqmcWalker;

    using LatticeBase = Lattice::LatticeBase;
    using ModelBase = Model::ModelBase;
    using MeasureHandler = Measure::MeasureHandler;


    // ----------------------- Interface class QuantumMonteCarlo::DqmcInitializer ------------------------
    class DqmcInitializer {
        public:

            static void initial_modules  ( LatticeBase& lattice, 
                                           ModelBase& model, 
                                           DqmcWalker& walker,
                                           MeasureHandler& meas_handler );

            static void initial_dqmc     ( LatticeBase& lattice, 
                                           ModelBase& model,
                                           DqmcWalker& walker,
                                           MeasureHandler& meas_handler );
    };

} // namespace QuantumMonteCarlo


#endif // DQMC_INITIALIZER_H