#include "hubbard.h"
#include "detQMC.h"
#include "options.h"

/**
 *  todo:
 *   1. get params from command lines (done)
 *   2. equal-time measurements of momentum distribution and spin-spin correlation (done)
 *   3. bin measurements (done)
 *   4. time-displaced green function and measurements (done)
 *   5. Stochastic Analytic Continuation (SAC) to obtain fermion spectrum function (missing)
 *   6. Check-board decomposition (missing)
 *   7. ******** Modify command console output ******** (done)
 *   8. attractive interaction U < 0 (done)
 *   9. ...
 */


/**  The Main Program */
int main(int argc, char* argv[]) {

    /* model and controlling params */
    int ll = 4;
    int lt = 80;
    double beta = 4.0;
    double t = 1.0;
    double u = 4.0;
    double mu = 0.0;

    int nwrap = 10;
    int nwarm = ceil(4 * ll * ll * beta);

    int nbin = 20;
    int nsweep = 100;
    int nBetweenBins = 10;

    std::string filename_eqtime = "../results/meas-eqtime.txt";
    std::string filename_dynamic = "../results/meas-dynamic.txt";
    bool bool_append = true;
    bool bool_display_process = false;

    bool bool_warm_up = true;
    bool bool_measure_eqtime = true;
    bool bool_measure_dynamic = true;
    
    /* read params from command line */
    getMyArgs(argc, argv, ll, lt, beta, t, u, mu, nwrap, nwarm, nbin, nsweep, nBetweenBins, filename_eqtime, filename_dynamic, bool_append, bool_measure_eqtime, bool_measure_dynamic);

    /* dqmc simulation */
    detQMC dqmc;

    for (double U = 0.0; U <= 3.0; U += 1.0) {
        bool_append = (U != 0);

        dqmc.set_Model_Params(ll, lt, beta, t, -U, mu, nwrap);

        dqmc.set_MC_Params(nwarm, nbin, nsweep, nBetweenBins);

        dqmc.set_bool_Params(bool_warm_up, bool_measure_eqtime, bool_measure_dynamic);

        dqmc.set_Momentum_q(M_PI, M_PI);

        dqmc.initialMeasure();

        dqmc.runQMC(bool_display_process);

        dqmc.analyseStats();

        dqmc.printStats();

        if (bool_measure_eqtime)
            dqmc.output_Stats_eqtime(filename_eqtime, bool_append);

        if (bool_measure_dynamic)
            dqmc.output_Stats_dynamic(filename_dynamic, bool_append);
    }

    /*
    dqmc.set_MC_Params(nwarm, nbin, nsweep, nBetweenBins);
    dqmc.set_Model_Params(ll, lt, beta, t, u, mu, nwrap);

    for (int i = 0; i <= ll/2; ++i) {
        // crystal momentum (qx, qy)
        const double qx = 2 * M_PI / ll * i;
        const double qy = 2 * M_PI / ll * i;

        dqmc.set_Momentum_q(qx, qy);

        bool_warm_up = (i == 0);
        dqmc.set_bool_Params(bool_warm_up, bool_measure_eqtime, bool_measure_dynamic);

        dqmc.initialMeasure();

        dqmc.runQMC(bool_display_process);

        dqmc.analyseStats();

        dqmc.printStats();

        bool_append = (i != 0);
        if (bool_measure_eqtime)
            dqmc.output_Stats_eqtime(filename_eqtime, bool_append);
    }
    */
}
