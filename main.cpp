#include "hubbard.h"
#include "detQMC.h"

/**
 *  todo:
 *   1. get params from command lines (boost) (missing)
 *   2. measurements of structure factor (missing)
 *   3. real-frequency information of electrons (missing)
 *   4. ...
 */

/**  The Main Program */
int main(int argc, char* argv[]) {

    int ll = 4;
    int lt = 80;
    double beta = 4.0;
    double t = 1.0;
    double U = 4.0;
    double mu = 0.0;

    int nwrap = 10;
    int nwarm = 4*ll*ll*beta;
    int nsweep = 200;

    std::string filename = "TestAndTest.dat";

    detQMC dqmc;

    dqmc.set_Model_Params(ll, lt, beta, t, U, mu, nwrap);

    dqmc.set_MC_Params(nwarm, nsweep);

    dqmc.runQMC(false);

    dqmc.printStats();

    dqmc.outputStats(filename, true);

}
