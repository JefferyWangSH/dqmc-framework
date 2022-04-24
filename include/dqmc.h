#ifndef DQMC_H
#define DQMC_H
#pragma once


/**
  *  
  *   
  */


#include <chrono>


namespace Model { class ModelBase; }
namespace Lattice { class LatticeBase; }
namespace Measure { class MeasureHandler; }


namespace QuantumMonteCarlo {

    // forward declaration
    class DqmcWalker;

    using ModelBase = Model::ModelBase;
    using LatticeBase = Lattice::LatticeBase;
    using MeasureHandler = Measure::MeasureHandler;
    

    // ------------------------ Pure interface class QuantumMonteCarlo::Dqmc -----------------------
    class Dqmc {
        
        public:

            static bool m_show_progress_bar;
            static void show_progress_bar( bool show_progress_bar );

            static std::chrono::steady_clock::time_point m_begin_time, m_end_time;

            static const double timer();
            
            static void thermalize           ( DqmcWalker& walker, 
                                               ModelBase& model,
                                               LatticeBase& lattice,  
                                               MeasureHandler& meas_handler );
            
            static void measure              ( DqmcWalker& walker, 
                                               ModelBase& model,
                                               LatticeBase& lattice,  
                                               MeasureHandler& meas_handler );

            static void analyse              ( MeasureHandler& meas_handler );


        private:

            static void sweep_forth_and_back ( DqmcWalker& walker, 
                                               ModelBase& model,
                                               LatticeBase& lattice, 
                                               MeasureHandler& meas_handler );

    };

}

#endif // DQMC_H
