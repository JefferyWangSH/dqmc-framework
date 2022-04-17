#include "dqmc_initializer.h"
#include "dqmc_walker.h"
#include "lattice/lattice_base.h"
#include "model/model_base.h"
#include "measure/measure_handler.h"


namespace QuantumMonteCarlo {

    void DqmcInitializer::initial_modules(  LatticeBase& lattice, 
                                            ModelBase& model, 
                                            DqmcWalker& walker,
                                            MeasureHandler& meas_handler )
    {
        // make sure that the params are setup correctly in advance,
        // and the orders of initializations below are important.

        // initialize lattice module
        lattice.initial();

        // initialize MeasureHandler module
        meas_handler.initial();

        // initialize dqmcWalker module
        walker.initial( lattice, meas_handler );

        // initialize model module
        model.initial( lattice, walker );
    }


    void DqmcInitializer::initial_dqmc( LatticeBase& lattice, 
                                        ModelBase& model, 
                                        DqmcWalker& walker,
                                        MeasureHandler& meas_handler )
    {
        // this subroutine should be called after the initial 
        // configuration of the bosonic fields have been determined, 
        // either randomly initialized or read from a input config file.
        // SvdStack class are initialized and the greens functions 
        // for the initial bosonic fields are computed in this function.
        walker.initial_svd_stacks( lattice, model );
        walker.initial_greens_function();
        walker.initial_config_sign();
    }


} // namespace QuantumMonteCarlo