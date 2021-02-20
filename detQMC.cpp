#include "detQMC.h"

void detQMC::set_Model_Params(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap) {
    // half-filled case only
    assert(mu == 0);
    Hubbard newHubb(ll, lt, beta, t, Uint, mu, nwrap);
    hubb = newHubb;
    this->nwrap = nwrap;
}

void detQMC::set_MC_Params(int nwarm, int nsweep) {
    this->nwarm = nwarm;
    this->nsweep = nsweep;
}

void detQMC::clearStats() {
    nn = 0;
    DoubleOccu = 0.0;
    KineticEnergy = 0.0;
    StructFactor = 0.0;
}

void detQMC::runQMC(bool bool_display_process) {

    clearStats();

    // record current time
    begin_t = clock();

    // thermalization process
    for (int nwm = 0; nwm < nwarm/2; ++nwm) {
        if ( nwm % 5 == 0 && bool_display_process)
        {
            std::cout << "warm up sweep: "<< 2 * nwm << std::endl;
        }
        sweep_BackAndForth(false);
    }
    // measuring process
    for (int nsw = 0; nsw < nsweep/2; ++nsw) {
        if ( nsw % 5 == 0 && bool_display_process)
        {
            std::cout << "measuring sweep: "<< 2 * nsw << std::endl;
        }
        sweep_BackAndForth(true);
    }

    end_t = clock();
}

void detQMC::sweep_BackAndForth(bool doMeasure) {
    // sweep back and forth
    hubb.sweep_0_to_beta(nwrap);
    if (doMeasure) { measure();}
    hubb.sweep_beta_to_0(nwrap);
    if (doMeasure) { measure();}
}

void detQMC::measure() {
    // read model params from hubbard
    int nx = hubb.ll, ny = hubb.ll;
    int lt = hubb.lt;
    double t = hubb.t;

    for (int l = 0; l < lt; ++l) {
        // loop for all lattice sites
        for (int x = 0; x < nx; ++x) {
            for (int y = 0; y < ny; ++y) {
                nn += 1;

                KineticEnergy += 2*t * hubb.vecGreenU[l](x + nx * y,((x+1) % nx) + nx * y)
                        + 2*t * hubb.vecGreenU[l](x + nx * y,x + nx * ((y+1) % ny));
                KineticEnergy += 2*t * hubb.vecGreenD[l](x + nx * y,((x+1) % nx) + nx * y)
                        + 2*t * hubb.vecGreenD[l](x + nx * y,x + nx * ((y+1) % ny));
                DoubleOccu += (1 - hubb.vecGreenU[l](x + nx * y, x + nx * y))
                        *(1 - hubb.vecGreenD[l](x + nx * y, x + nx * y));
                // todo: StructFactor
            }
        }
    }
}

void detQMC::printStats() {

    double time = (double)(end_t - begin_t)/CLOCKS_PER_SEC;
    int minute = floor(time / 60);
    double sec = time - 60 * minute;

    std::cout << "ll:  " << hubb.ll << std::endl
              << "lt:  " << hubb.lt << std::endl
              << "beta: " << hubb.beta << std::endl
              << "U/t:  " << hubb.Uint / hubb.t << std::endl
              << "mu:   " << hubb.mu << std::endl
              << "nwrap:  " << nwrap << std::endl
              << "Order param:    " << DoubleOccu / nn << std::endl
              << "Kinetic energy: " << KineticEnergy / nn << std::endl
              << "time cost:      " << minute << " min "
              << std::setiosflags(std::ios::fixed) << std::setprecision(1) << sec << " s" << std::endl
              << std::endl;
    std::cout.unsetf(std::ios::adjustfield|std::ios::basefield|std::ios::floatfield);
}

void detQMC::outputStats(const std::string &filename, bool bool_Append) {

    std::ofstream outfile;
    if (bool_Append) {
        outfile.open(filename, std::ios::out | std::ios::app);
    }
    else {
        outfile.open(filename, std::ios::out | std::ios::trunc);
    }
    outfile << std::setiosflags(std::ios::right)
            << std::setw(15) << hubb.Uint / hubb.t
            << std::setw(15) << DoubleOccu / nn
            << std::setw(15) << KineticEnergy / nn
            << std::endl;
    outfile.close();
}


