#include "detQMC.h"
#include "ProgressBar.hpp"
#include <cmath>

void detQMC::set_Model_Params(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap) {
    // half-filled case only, free of sign problem
    assert(mu == 0);
    Hubbard newHubb(ll, lt, beta, t, Uint, mu, nwrap);
    hubb = newHubb;
    this->nwrap = nwrap;
}

void detQMC::set_MC_Params(int nwarm, int nbin, int nsweep, int nBetweenBins) {
    this->nwarm = nwarm;
    this->nbin = nbin;
    this->nsweep = nsweep;
    this->nBetweenBins = nBetweenBins;
}

void detQMC::set_bool_Params(bool bool_warm_up, bool bool_measure_eqtime, bool bool_measure_dynamic) {
    this->bool_warm_up = bool_warm_up;
    this->bool_measure_eqtime = bool_measure_eqtime;
    this->bool_measure_dynamic = bool_measure_dynamic;
}

void detQMC::set_Momentum_q(double qx, double qy) {
    q = (vecXd(2) << qx, qy).finished();
}

void detQMC::printParams() {
    std::cout << "===========================================================================" << std::endl;
    std::cout << "  Simulation Parameters: " << std::endl
              << "    ll:  " << hubb.ll << std::endl
              << "    lt:  " << hubb.lt << std::endl
              << "    beta: " << hubb.beta << std::endl
              << "    U/t:  " << hubb.Uint / hubb.t << std::endl
              << "    mu:   " << hubb.mu << std::endl
              << "    q:    " << "(" << q(0) << ", "<< q(1) << ")" << std::endl
              << "    nwrap:  " << nwrap << std::endl;
    std::cout << "===========================================================================" << std::endl;
}

void detQMC::initialMeasure() {
    // initialize bins of observables
    if (bool_measure_eqtime) {
        obs_bin_eqtime["DoubleOccu"].reserve(nbin);
        obs_bin_eqtime["KineticEnergy"].reserve(nbin);
        obs_bin_eqtime["StructFactor"].reserve(nbin);
        obs_bin_eqtime["MomentumDist"].reserve(nbin);
        obs_bin_eqtime["localSpinCorr"].reserve(nbin);
    }

    if (bool_measure_dynamic) {
        for (int bin = 0; bin < nbin; ++bin) {
            for (int l = 0; l < hubb.lt; ++l) {
                obs_bin_gt0[bin][l] = matXd::Zero(hubb.ls, hubb.ls);
            }
        }

        for (int l = 0; l < hubb.lt; ++l) {
            obs_mean_gt0_k[l] = 0.0;
            obs_err_gt0_k[l] = 0.0;
            vec_gt0_tau[l] = matXd::Identity(hubb.ls, hubb.ls);
        }
    }
}

void detQMC::clearStats() {
    if (bool_measure_eqtime) {
        n_equal_time = 0;
        DoubleOccu = 0.0;
        KineticEnergy = 0.0;
        StructFactor = 0.0;
        MomentumDist = 0.0;
        localSpinCorr = 0.0;
    }

    if (bool_measure_dynamic) {
        n_time_displaced = 0;
        for (int l = 0; l < hubb.lt; ++l) {
            vec_gt0_tau[l] = matXd::Zero(hubb.ls, hubb.ls);
        }
    }
}

void detQMC::runQMC(bool bool_display_process) {

    clearStats();

    // record current time
    begin_t = clock();

    if (bool_warm_up) {
        // thermalization process

        // progress bar
        progresscpp::ProgressBar progressBar(nwarm/2, 40, '#', '-');

        for (int nwm = 1; nwm <= nwarm/2; ++nwm) {
            sweep_BackAndForth(false, false);
            ++progressBar;

            if ( nwm % 10 == 0 && bool_display_process) {
                std::cout << "Warm-up progress:   ";
                progressBar.display();
            }
        }

        if (bool_display_process) {
            std::cout << "Warm-up progress:   ";
            progressBar.done();
        }
    }

    if (bool_measure_eqtime || bool_measure_dynamic) {
        // measuring process

        progresscpp::ProgressBar progressBar(nbin * nsweep / 2, 40, '#', '-');

        for (int bin = 1; bin <= nbin; ++bin) {
            for (int nsw = 1; nsw <= nsweep/2; ++nsw) {
                sweep_BackAndForth(bool_measure_eqtime, bool_measure_dynamic);
                ++progressBar;

                if ( nsw % 10 == 0 && bool_display_process) {
                    std::cout << "Measuring progress: ";
                    progressBar.display();
                }
            }

            // analyse statistical data
            normalizeStats();

            write_Stats_to_bins(bin);

            clearStats();

            // avoid correlation between bins
            for (int n_bw = 0; n_bw < nBetweenBins; ++n_bw) {
                sweep_BackAndForth(false, false);
            }
        }

        if (bool_display_process) {
            std::cout << "Measuring progress: ";
            progressBar.done();
        }
    }

    end_t = clock();
}

void detQMC::sweep_BackAndForth(bool bool_eqtime, bool bool_dynamic) {

    // sweep forth from 0 to beta
    if (!bool_dynamic) {
        hubb.sweep_0_to_beta(nwrap);
    }
    else {
        hubb.sweep_0_to_beta_displaced(nwrap);
        measure_time_displaced();
    }
    if (bool_eqtime) {
        measure_equal_time();
    }

    // sweep back from beta to 0
    hubb.sweep_beta_to_0(nwrap);
    if (bool_eqtime) {
        measure_equal_time();
    }
}

void detQMC::measure_equal_time() {

    for (int l = 0; l < hubb.lt; ++l) {
        matXd gu = hubb.vecGreenU[l];
        matXd gd = hubb.vecGreenD[l];

        n_equal_time++;
        meas_Double_Occu(gu, gd);
        meas_Kinetic_Energy(gu, gd);
        meas_Struct_Factor(gu, gd, q);
        meas_Momentum_Dist(gu, gd, q);
        meas_local_Spin_Corr(gu, gd);
    }
}

void detQMC::measure_time_displaced() {
    n_time_displaced++;
    for (int l = 0; l < hubb.lt; ++l) {
        /* factor 2 comes from two spin states */
        vec_gt0_tau[l] += (hubb.vecGreen_t0_up[l] + hubb.vecGreen_t0_dn[l]) / 2;
    }
}

void detQMC::normalizeStats() {

    if (bool_measure_eqtime) {
        DoubleOccu /= hubb.ls * n_equal_time;
        KineticEnergy /= hubb.ls * n_equal_time;
        StructFactor /= hubb.ls * hubb.ls * n_equal_time;
        MomentumDist /= n_equal_time;
        localSpinCorr /= n_equal_time;
    }

    if (bool_measure_dynamic) {
        for (int l = 0; l < hubb.lt; ++l) {
            vec_gt0_tau[l] /= n_time_displaced;
        }
    }
}

void detQMC::write_Stats_to_bins(int bin) {
    if (bool_measure_eqtime) {
        obs_bin_eqtime["DoubleOccu"][bin] = DoubleOccu;
        obs_bin_eqtime["KineticEnergy"][bin] = KineticEnergy;
        obs_bin_eqtime["StructFactor"][bin] = StructFactor;
        obs_bin_eqtime["MomentumDist"][bin] = MomentumDist;
        obs_bin_eqtime["localSpinCorr"][bin] = localSpinCorr;
    }

    if (bool_measure_dynamic) {
        for (int l = 0; l < hubb.lt; ++l) {
            obs_bin_gt0[bin][l] = vec_gt0_tau[l];
        }
    }
}

void detQMC::analyse_equal_time_Stats(const std::string& obs) {
    assert(obs_bin_eqtime.count(obs) == 1);

    // clear data of previous statistics
    obs_mean_eqtime[obs] = 0;
    obs_err_eqtime[obs] = 0;

    for (int i = 0; i < nbin; ++i) {
        obs_mean_eqtime[obs] += obs_bin_eqtime[obs][i];
        obs_err_eqtime[obs] += pow(obs_bin_eqtime[obs][i], 2);
    }

    obs_mean_eqtime[obs] /= nbin;
    obs_err_eqtime[obs] /= nbin;
    obs_err_eqtime[obs] = pow(obs_err_eqtime[obs]-pow(obs_mean_eqtime[obs], 2), 0.5) / pow(nbin - 1, 0.5);
}

void detQMC::analyse_timeDisplaced_Stats(const vecXd &p) {

    // calculate the Matsubara green function in momentum space and record data.
    for (int l = 0; l < hubb.lt; ++l) {
        // clear previous statistics
        obs_mean_gt0_k[l] = 0.0;
        obs_err_gt0_k[l] = 0.0;

        for (int bin = 0; bin < nbin; ++bin) {
            const matXd Matsubara = obs_bin_gt0[bin][l];
            double tmpFourier = 0.0;

            for (int xi = 0; xi < hubb.ll; ++xi) {
                for (int yi = 0; yi < hubb.ll; ++yi) {
                    for (int xj = 0; xj < hubb.ll; ++xj) {
                        for (int yj = 0; yj < hubb.ll; ++yj) {
                            const int i = xi + hubb.ll * yi;
                            const int j = xj + hubb.ll * yj;
                            const vecXd r = (vecXd(2) << (xi-xj), (yi-yj)).finished();
                            // TODO: check complex or real
                            tmpFourier += cos(-r.dot(p)) * Matsubara(j, i) / hubb.ls;
                        }
                    }
                }
            }

            obs_mean_gt0_k[l] += tmpFourier;
            obs_err_gt0_k[l] += tmpFourier * tmpFourier;
        }
        obs_mean_gt0_k[l] /= nbin;
        obs_err_gt0_k[l] /= nbin;
        obs_err_gt0_k[l] = pow(obs_err_gt0_k[l]-pow(obs_mean_gt0_k[l], 2), 0.5) / pow(nbin - 1, 0.5);
    }
}

void detQMC::analyseStats() {
    if (bool_measure_eqtime) {
        analyse_equal_time_Stats("DoubleOccu");
        analyse_equal_time_Stats("KineticEnergy");
        analyse_equal_time_Stats("StructFactor");
        analyse_equal_time_Stats("MomentumDist");
        analyse_equal_time_Stats("localSpinCorr");
    }

    if (bool_measure_dynamic) {
        // FIXME: support sequence of momentum q
        analyse_timeDisplaced_Stats(q);
    }
}

void detQMC::printStats() {

    double time = (double)(end_t - begin_t)/CLOCKS_PER_SEC;
    int minute = floor(time / 60);
    double sec = time - 60 * minute;

    if (bool_measure_eqtime) {
        std::cout.precision(8);
        std::cout << std::endl;
        std::cout << "  Equal-time Measurements: " << std::endl
                  << "    Double occu:      " << obs_mean_eqtime["DoubleOccu"] << "     err: " << obs_err_eqtime["DoubleOccu"] << std::endl
                  << "    Kinetic energy:   " << obs_mean_eqtime["KineticEnergy"] << "     err: " << obs_err_eqtime["KineticEnergy"] << std::endl
                  << "    Momentum dist:    " << obs_mean_eqtime["MomentumDist"] << "     err: " << obs_err_eqtime["MomentumDist"] << std::endl
                  << "    local Spin corr:  " << obs_mean_eqtime["localSpinCorr"] << "     err: " << obs_err_eqtime["localSpinCorr"] << std::endl
                  << "    Structure Factor: " << obs_mean_eqtime["StructFactor"] << "     err: " << obs_err_eqtime["StructFactor"] << std::endl
                  << std::endl;
        std::cout.precision(-1);
    }

    if (bool_measure_dynamic) {
        std::cout << "  Time-displaced Measurements: " << std::endl
                  << "    Matsubara Green's function in momentum space:  see in file" << std::endl
                  << std::endl;
    }

    std::cout << "  Time Cost:      " << minute << " min " << sec << " s" << std::endl;

    std::cout << "===========================================================================" << std::endl
              << std::endl;
}

void detQMC::output_Stats_eqtime(const std::string &filename, bool bool_Append) {
    assert(bool_measure_eqtime);

    std::ofstream outfile;
    if (bool_Append) {
        outfile.open(filename, std::ios::out | std::ios::app);
    }
    else {
        outfile.open(filename, std::ios::out | std::ios::trunc);
    }
    outfile << std::setiosflags(std::ios::right)
            << std::setw(15) << hubb.Uint / hubb.t
            << std::setw(15) << hubb.beta
            << std::setw(15) << obs_mean_eqtime["DoubleOccu"]
            << std::setw(15) << obs_mean_eqtime["KineticEnergy"]
            << std::setw(15) << obs_mean_eqtime["StructFactor"]
            << std::setw(15) << obs_mean_eqtime["MomentumDist"]
            << std::setw(15) << obs_mean_eqtime["localSpinCorr"]
            << std::setw(15) << obs_err_eqtime["DoubleOccu"]
            << std::setw(15) << obs_err_eqtime["KineticEnergy"]
            << std::setw(15) << obs_err_eqtime["StructFactor"]
            << std::setw(15) << obs_err_eqtime["MomentumDist"]
            << std::setw(15) << obs_err_eqtime["localSpinCorr"]
            << std::setw(15) << q(0)
            << std::setw(15) << q(1)
            << std::endl;
    outfile.close();
}

void detQMC::output_Stats_dynamic(const std::string& filename, bool bool_Append) {
    assert(bool_measure_dynamic);

    std::ofstream outfile;
    if (bool_Append) {
        outfile.open(filename, std::ios::out | std::ios::app);
    }
    else {
        outfile.open(filename, std::ios::out | std::ios::trunc);
    }

    outfile << std::setiosflags(std::ios::right)
            << "Momentum k: "<< "(" << q(0) << ", "<< q(1) << ")" << std::endl;

    for (int l = 1; l <= hubb.lt; ++l) {
        outfile << std::setw(15) << l
                << std::setw(15) << obs_mean_gt0_k[l-1]
                << std::setw(15) << obs_err_gt0_k[l-1]
                << std::setw(15) << obs_err_gt0_k[l-1] / obs_mean_gt0_k[l-1]
                << std::endl;
    }
    outfile.close();
}

void detQMC::meas_Double_Occu(const matXd& gu, const matXd& gd) {
    const int ls = hubb.ls;

    for (int i = 0; i < ls; ++i) {
        const double doubleoccu = (1 - gu(i,i)) * (1 - gd(i,i));
        DoubleOccu += doubleoccu;
    }
}

void detQMC::meas_Kinetic_Energy(const matXd& gu, const matXd& gd) {
    const int ll = hubb.ll;
    const double t = hubb.t;

    for (int x = 0; x < ll; ++x) {
        for (int y = 0; y < ll; ++y) {
            const double kinetic = 2 * t * (gu(x + ll*y, ((x+1)%ll) + ll*y) + gu(x + ll*y, x + ll*((y+1)%ll)))
                    + 2 * t * (gd(x + ll*y, ((x+1)%ll) + ll*y) + gd(x + ll*y, x + ll*((y+1)%ll)));
            KineticEnergy += kinetic;
        }
    }
}

void detQMC::meas_Momentum_Dist(const matXd &gu, const matXd &gd, const vecXd& p) {
    const int ll = hubb.ll;
    double tmpfourier = 0.0;

    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            for (int xj = 0; xj < ll; ++xj) {
                for (int yj = 0; yj < ll; ++yj) {
                    const int i = xi + ll * yi;
                    const int j = xj + ll * yj;
                    const vecXd r = (vecXd(2) << (xi-xj), (yi-yj)).finished();
                    tmpfourier += cos(-r.dot(p)) * (gu(j, i) + gd(j, i));
                }
            }
        }
    }
    MomentumDist += 1 - tmpfourier/2/hubb.ls;
}

void detQMC::meas_local_Spin_Corr(const matXd &gu, const matXd &gd) {
    const int ls = hubb.ls;
    double  onsitecorr = 0.0;

    for (int i = 0; i < ls; ++i) {
        onsitecorr += gu(i,i) + gd(i,i) - 2*gu(i,i)*gd(i,i);
    }
    localSpinCorr += onsitecorr/ls;
}


void detQMC::meas_Struct_Factor(const matXd &gu, const matXd &gd, const vecXd &p) {
    const int ll = hubb.ll;
    const int ls = hubb.ls;

    /**  gu(i,j) = < c_i c^+_j >
     *  guc(i,j) = < c^+_i c_j > */
    matXd guc = matXd::Identity(ls, ls);
    matXd gdc = matXd::Identity(ls, ls);

    // get guc and gdc
    for (int i = 0; i < ls; ++i) {
        for (int j = 0; j < ls; ++j) {
            guc(j,i) = - gu(i,j);
            gdc(j,i) = - gd(i,j);
        }
        guc(i,i)++;
        gdc(i,i)++;
    }

    // loop for site i, j
    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            for (int xj = 0; xj < ll; ++xj) {
                for (int yj = 0; yj < ll; ++yj) {
                    const int i = xi + ll * yi;
                    const int j = xj + ll * yj;
                    const vecXd r = (vecXd(2) << (xi-xj), (yi-yj)).finished();
                    const double factor = cos(-r.dot(p));
                    /** factor 4 comes from spin 1/2 */
                    const double structfactor = factor / 4 * (guc(i,i)*guc(j,j) + guc(i,j)*gu(i,j)
                                                            + gdc(i,i)*gdc(j,j) + gdc(i,j)*gd(i,j)
                                                            - gdc(i,i)*guc(j,j) - guc(i,i)*gdc(j,j));
                    StructFactor += structfactor;
                }
            }
        }
    }
}

detQMC::~detQMC() {
    std::cout << std::endl << "The simulation was done :)" << std::endl;
}
