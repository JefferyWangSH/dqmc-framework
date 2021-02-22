#include "hubbard.h"
#include "detQMC.h"
#include "options.h"


/**
 *  todo:
 *   1. get params from command lines (boost) (missing, temporarily using options.h)
 *   2. measurements of structure factor (missing)
 *   3. real-frequency information of electrons (missing)
 *   4. ...
 */

/**  The Main Program */
int main(int argc, char* argv[]) {

    /* model and controlling params */
    int ll = 4;
    int lt = 80;
    double beta = 4.0;
    double t = 1.0;
    double U = 4.0;
    double mu = 0.0;

    int nwrap = 10;
    int nwarm = 4 * ll * ll * beta;
    int nsweep = 200;

    std::string filename = "output.txt";
    bool bool_append = true;

    /* read params from command line */
    getMyArgs(argc, argv, ll, lt, beta, t, U, mu, nwrap, nwarm, nsweep, filename, bool_append);

    /* dqmc simulation */
    detQMC dqmc;

    dqmc.set_Model_Params(ll, lt, beta, t, U, mu, nwrap);
    dqmc.set_MC_Params(nwarm, nsweep);

    dqmc.runQMC(true);

    dqmc.printStats();
    dqmc.outputStats(filename, bool_append);

}
