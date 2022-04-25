#ifndef DQMC_H
#define DQMC_H
#pragma once


/**
  *  This header file defines QuantumMonteCarlo::Dqmc class for the organizations of the dqmc program.
  *  Top-level designs for the dqmc simulation, e.g. thermalization, measurements and data analysis,
  *  are implemented as static member functions of this Dqmc class.
  *  Some other useful tools, such as timer and progress bar, 
  *  are also provided under the namespace QuantumMonteCarlo::Dqmc .
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
    

    // -------------------------------- Pure interface class QuantumMonteCarlo::Dqmc --------------------------------
    class Dqmc {
        
        public:

            // ---------------------------------------- Useful tools ------------------------------------------
            
            // set up whether to show the process bar or not
            static void show_progress_bar( bool show_progress_bar );

            // set up the format of the progress bar
            static void progress_bar_format( unsigned int width, char complete, char incomplete );
            
            // return the duration time of the dqmc process, e.g thermalization or measurements
            static const double timer();

            
            // ------------------------------------ Crucial Dqmc routines -------------------------------------
            
            // thermalization of the field configurations
            static void thermalize           ( DqmcWalker& walker, 
                                               ModelBase& model,
                                               LatticeBase& lattice,  
                                               MeasureHandler& meas_handler );
            
            // Monte Carlo updates and measurments
            static void measure              ( DqmcWalker& walker, 
                                               ModelBase& model,
                                               LatticeBase& lattice,  
                                               MeasureHandler& meas_handler );

            // analyse the measured data
            static void analyse              ( MeasureHandler& meas_handler );


        private:

            // ------------------------------------ Private subroutines -------------------------------------
            // declarations of static members
            static bool m_show_progress_bar;
            static unsigned int m_progress_bar_width;
            static char m_progress_bar_complete_char, m_progress_bar_incomplete_char;
            static std::chrono::steady_clock::time_point m_begin_time, m_end_time;

            // start the timer
            static void timer_begin();
            
            // end the timer
            static void timer_end();

            // sweep and update the field configurations 
            // from 0 to beta and back from beta to 0
            // do the measurements if needed
            static void sweep_forth_and_back ( DqmcWalker& walker, 
                                               ModelBase& model,
                                               LatticeBase& lattice, 
                                               MeasureHandler& meas_handler );

    };

}

#endif // DQMC_H
