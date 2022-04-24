#include "dqmc.h"
#include "dqmc_walker.h"
#include "model/model_base.h"
#include "lattice/lattice_base.h"
#include "measure/measure_handler.h"


namespace QuantumMonteCarlo {

    void Dqmc::sweep_forth_and_back( DqmcWalker& walker, 
                                     ModelBase& model, 
                                     LatticeBase& lattice, 
                                     MeasureHandler& meas_handler )
    {
        // sweep forth from 0 to beta
        if ( meas_handler.isDynamic() ) {
            walker.sweep_for_dynamic_greens(model);
            meas_handler.dynamic_measure(model, lattice);
        }
        else {
            walker.sweep_from_0_to_beta(model);
            if ( meas_handler.isEqualTime() ) {
                meas_handler.equaltime_measure(model, lattice);
            }
        }

        // sweep back from beta to 0
        walker.sweep_from_beta_to_0(model);
        if ( meas_handler.isEqualTime() ) {
            meas_handler.equaltime_measure(model, lattice);
        }
    }


    void Dqmc::thermalize( DqmcWalker& walker, 
                           ModelBase& model,
                           LatticeBase& lattice,  
                           MeasureHandler& meas_handler ) 
    {
        if ( meas_handler.isWarmUp() ) {
            for ( auto sweep = 1; sweep <= meas_handler.WarmUpSweeps()/2; sweep++ ) {
                // sweep forth and back without measuring
                walker.sweep_from_0_to_beta(model);
                walker.sweep_from_beta_to_0(model);
            }
        }
    }

    
    void Dqmc::measure( DqmcWalker& walker, 
                        ModelBase& model,
                        LatticeBase& lattice,  
                        MeasureHandler& meas_handler ) 
    {   
        if ( meas_handler.isEqualTime() || meas_handler.isDynamic() ) {
            for ( auto bin = 0; bin < meas_handler.BinsNum(); ++bin ) {
                for ( auto sweep = 1; sweep <= meas_handler.BinsSize()/2; ++sweep ) {
                    // update and measure
                    sweep_forth_and_back(walker, model, lattice, meas_handler);
                }

                // store the collected data in the MeasureHandler
                meas_handler.normalize_stats();
                meas_handler.write_stats_to_bins(bin);
                meas_handler.clear_temporary();

                // avoid correlations between adjoining bins
                for ( auto sweep = 0; sweep < meas_handler.SweepsBetweenBins()/2; ++sweep ) {
                    walker.sweep_from_0_to_beta(model);
                    walker.sweep_from_beta_to_0(model);
                }
            }
        }
    }


    void Dqmc::analyse( MeasureHandler& meas_handler )
    {
        // analyse the collected data after the measuring process
        meas_handler.analyse_stats();
    }
            

} // namespace QuantumMonteCarlo
