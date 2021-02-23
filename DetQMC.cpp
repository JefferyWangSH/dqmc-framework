#include "detQMC.h"
#include <complex>
#include <cmath>

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

void detQMC::set_Momentum_q(double qx, double qy) {
    q = (vecXd(2) << qx, qy).finished();
}

void detQMC::clearStats() {
    nn = 0;
    DoubleOccu = 0.0;
    KineticEnergy = 0.0;
    StructFactor = 0.0;
    MomentumDist = 0.0;
    localSpinCorr = 0.0;
}

void detQMC::runQMC(bool bool_warm, bool bool_display_process) {

    clearStats();

    // record current time
    begin_t = clock();

    if (bool_warm)
    {
        // thermalization process
        for (int nwm = 0; nwm < nwarm/2; ++nwm) {
            if ( nwm % 5 == 0 && bool_display_process)
            {
                std::cout << "warm up sweep: "<< 2 * nwm << std::endl;
            }
            sweep_BackAndForth(false);
        }
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
    if (doMeasure) { measure_equal_time();}
    hubb.sweep_beta_to_0(nwrap);
    if (doMeasure) { measure_equal_time();}
}

void detQMC::measure_equal_time() {

    for (int l = 0; l < hubb.lt; ++l) {
        matXd gu = hubb.vecGreenU[l];
        matXd gd = hubb.vecGreenD[l];

        /* todo: update 'measurements' in dqmc note */
        nn++;
        DoubleOccu += meas_DoubleOccu(gu, gd);
        KineticEnergy += meas_KineticEnergy(gu, gd);
        MomentumDist += meas_MomentumDist(gu, gd, q);
        localSpinCorr += meas_localSpinCorr(gu, gd);
        StructFactor += meas_StructFactor(gu, gd ,q);
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
              << "q:    " << "(" << q(0) << ", "<< q(1) << ")" << std::endl
              << "nwrap:  " << nwrap << std::endl
              << "Order param:      " << DoubleOccu / nn << std::endl
              << "Kinetic energy:   " << KineticEnergy / nn << std::endl
              << "Momentum dist:    " << MomentumDist / nn << std::endl
              << "local Spin corr:  " << localSpinCorr / nn << std::endl
              << "Structure Factor: " << StructFactor / nn << std::endl
              << "time cost:      " << minute << " min " << sec << " s" << std::endl
              << std::endl;
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
            << std::setw(15) << MomentumDist / nn
            << std::setw(15) << localSpinCorr / nn
            << std::setw(15) << StructFactor / nn
            << std::setw(15) << q(0)
            << std::setw(15) << q(1)
            << std::endl;
    outfile.close();
}

double detQMC::meas_DoubleOccu(const matXd& gu, const matXd& gd) const {
    const int ls = hubb.ls;
    double doubleoccu = 0.0;

    for (int i = 0; i < ls; ++i) {
        doubleoccu += (1 - gu(i,i)) * (1 - gd(i,i));
    }
    return doubleoccu/ls;
}

double detQMC::meas_KineticEnergy(const matXd& gu, const matXd& gd) const {
    const int ll = hubb.ll;
    const double t = hubb.t;
    double kinetic = 0.0;

    for (int x = 0; x < ll; ++x) {
        for (int y = 0; y < ll; ++y) {
            kinetic += 2 * t * (gu(x + ll*y, ((x+1)%ll) + ll*y) + gu(x + ll*y, x + ll*((y+1)%ll)));
            kinetic += 2 * t * (gd(x + ll*y, ((x+1)%ll) + ll*y) + gd(x + ll*y, x + ll*((y+1)%ll)));
        }
    }
    return kinetic/ll/ll;
}

double detQMC::meas_MomentumDist(const matXd &gu, const matXd &gd, const vecXd& p) const {
    const int ll = hubb.ll;
    double tmpfourier = 0.0;

    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            for (int xj = 0; xj < ll; ++xj) {
                for (int yj = 0; yj < ll; ++yj) {
                    const int i = xi + ll * yi;
                    const int j = xj + ll * yj;

                    vecXd r = (vecXd(2) << (xi-xj), (yi-yj)).finished();
                    std::complex<double> phase(0, -r.dot(p));
                    const double factor = exp(phase).real();

                    tmpfourier += factor * (gu(j, i) + gd(j, i));
                }
            }
        }
    }
    return 1 - tmpfourier/2/hubb.ls;
}

double detQMC::meas_localSpinCorr(const matXd &gu, const matXd &gd) const {
    const int ls = hubb.ls;
    double onsitecorr = 0.0;

    for (int i = 0; i < ls; ++i) {
        onsitecorr += gu(i,i) + gd(i,i) - 2*gu(i,i)*gd(i,i);
    }
    return onsitecorr/ls;
}

double detQMC::meas_StructFactor(const matXd &gu, const matXd &gd, const vecXd &p) const {
    const int ll = hubb.ll;
    double structfactor = 0.0;

    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            for (int xj = 0; xj < ll; ++xj) {
                for (int yj = 0; yj < ll; ++yj) {
                    const int i = xi + ll * yi;
                    const int j = xj + ll * yj;

                    vecXd r = (vecXd(2) << (xi-xj), (yi-yj)).finished();
                    std::complex<double> phase(0, -r.dot(p));
                    const double factor = exp(phase).real();

                    if (i == j)
                        structfactor += factor * ( gu(i,i) - 2 * gu(i,i) * gd(i,i) + gd(i,i) );
                    else
                        structfactor += factor * ( (gu(i,i)-gd(i,i)) * (gu(j,j)-gd(j,j)) - gu(j,i)*gu(i,j) - gd(j,i)*gd(i,j) );
                }
            }
        }
    }
    return structfactor/4/hubb.ls/hubb.ls;
}
