#include "hubbard.h"
#include "detQMC.h"
#include "options.h"


/**
 *  todo:
 *   1. get params from command lines (done options.h)
 *   2. measurements of (real-frequency) momentum distribution and spin correlation (missing)
 *   3. displaced greens function and measurements (missing)
 *   4. ...
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
    int nwarm = 4 * ll * ll * beta;
    int nsweep = 200;

    std::string filename = "output.txt";
    bool bool_append = true;

    /* read params from command line */
    getMyArgs(argc, argv, ll, lt, beta, t, u, mu, nwrap, nwarm, nsweep, filename, bool_append);

    /* dqmc simulation */
    detQMC dqmc;

    for (double U = 0.0; U <= 8.0; U += 0.25) {
        if (U == 0) { bool_append = false; }
        else { bool_append = true; }

        dqmc.set_Model_Params(ll, lt, beta, t, U, mu, nwrap);
        dqmc.runQMC(true, false);

        dqmc.printStats();
        dqmc.outputStats(filename, bool_append);
    }

    /*
    dqmc.set_MC_Params(nwarm, nsweep);
    dqmc.set_Model_Params(ll, lt, beta, t, u, mu, nwrap);

    for (int i = 0; i <= 30; ++i) {
        const double qx = M_PI/30 * i;
        const double qy = M_PI/30 * i;

        dqmc.set_Momentum_q(qx, qy);

        bool bool_warm = (i == 0);
        dqmc.runQMC(bool_warm, true);

        dqmc.printStats();
        bool_append = (i != 0);
        dqmc.outputStats(filename, bool_append);
    }*/



}

