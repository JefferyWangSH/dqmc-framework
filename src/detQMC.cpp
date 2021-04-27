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
    this->nsweep = nsweep;
    this->nBetweenBins = nBetweenBins;
    this->nbin = nbin;

    eqtimeMeasure.resize(nbin);
    dynamicMeasure.resize(nbin);
}

void detQMC::set_bool_Params(bool bool_warm_up, bool bool_measure_eqtime, bool bool_measure_dynamic) {
    this->bool_warm_up = bool_warm_up;
    this->bool_measure_eqtime = bool_measure_eqtime;
    this->bool_measure_dynamic = bool_measure_dynamic;
}

void detQMC::set_Momentum_q(double qx, double qy) {
    q = (Eigen::VectorXd(2) << qx, qy).finished();
    eqtimeMeasure.q = (Eigen::VectorXd(2) << qx, qy).finished();
    dynamicMeasure.q = (Eigen::VectorXd(2) << qx, qy).finished();
}

void detQMC::printParams() {
    std::cout << std::endl;
    std::cout << "==============================================================================" << std::endl;
    std::cout << "  Simulation Parameters: " << std::endl
              << "    ll:  " << hubb.ll << std::endl
              << "    lt:  " << hubb.lt << std::endl
              << "    beta: " << hubb.beta << std::endl
              << "    U/t:  " << hubb.Uint / hubb.t << std::endl
              << "    mu:   " << hubb.mu << std::endl
              << "    q:    " << "(" << q(0) << ", "<< q(1) << ")" << std::endl
              << "    nwrap:  " << nwrap << std::endl;
    std::cout << "==============================================================================" << std::endl;
}

void detQMC::initialMeasure() {
    // initialize bins of observables
    if (bool_measure_eqtime) {
        eqtimeMeasure.initial();
    }

    if (bool_measure_dynamic) {
        dynamicMeasure.initial(hubb);
    }
}

void detQMC::runQMC(bool bool_display_process) {

    // clear data
    if (bool_measure_eqtime) {
        eqtimeMeasure.clear();
    }

    if (bool_measure_dynamic) {
        dynamicMeasure.clear(hubb);
    }

    // record current time
    begin_t = std::chrono::steady_clock::now();

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
            if (bool_measure_eqtime) {
                eqtimeMeasure.normalizeStats(hubb);
                eqtimeMeasure.write_Stats_to_bins(bin);
                eqtimeMeasure.clear();
            }

            if (bool_measure_dynamic) {
                dynamicMeasure.normalizeStats(hubb);
                dynamicMeasure.write_Stats_to_bins(bin, hubb);
                dynamicMeasure.clear(hubb);
            }

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

    end_t = std::chrono::steady_clock::now();
}

void detQMC::sweep_BackAndForth(bool bool_eqtime, bool bool_dynamic) {

    // sweep forth from 0 to beta
    if (!bool_dynamic) {
        hubb.sweep_0_to_beta(nwrap);
    }
    else {
        hubb.sweep_0_to_beta_displaced(nwrap);
        dynamicMeasure.measure_time_displaced(hubb);
    }
    if (bool_eqtime) {
        eqtimeMeasure.measure_equal_time(hubb);
    }

    // sweep back from beta to 0
    hubb.sweep_beta_to_0(nwrap);
    // todo: hubb.sweep_beta_to_0_displaced
    if (bool_eqtime) {
        eqtimeMeasure.measure_equal_time(hubb);
    }
}

void detQMC::analyseStats() {
    if (bool_measure_eqtime) {
        eqtimeMeasure.analyseStats();
    }

    if (bool_measure_dynamic) {
        dynamicMeasure.analyse_timeDisplaced_Stats(hubb);
    }
}

void detQMC::printStats() {

    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count();

    const int minute = std::floor((double)time / 1000 / 60);
    const double sec = (double)time / 1000 - 60 * minute;

    if (bool_measure_eqtime) {
        std::cout.precision(8);
        std::cout << std::endl;
        std::cout << "  Equal-time Measurements: " << std::endl
                  << "    Double occu:      " << eqtimeMeasure.obs_mean_eqtime["DoubleOccu"]
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["DoubleOccu"] << std::endl
                  << "    Kinetic energy:   " << eqtimeMeasure.obs_mean_eqtime["KineticEnergy"]
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["KineticEnergy"] << std::endl
                  << "    Momentum dist:    " << eqtimeMeasure.obs_mean_eqtime["MomentumDist"]
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["MomentumDist"] << std::endl
                  << "    local Spin corr:  " << eqtimeMeasure.obs_mean_eqtime["localSpinCorr"]
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["localSpinCorr"] << std::endl
                  << "    Structure Factor: " << eqtimeMeasure.obs_mean_eqtime["StructFactor"]
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["StructFactor"] << std::endl;
        std::cout.precision(-1);
    }

    if (bool_measure_dynamic) {
        std::cout.precision(8);
        std::cout << std::endl;
        std::cout << "  Time-displaced Measurements: " << std::endl
                  << "    Dynamical correlation in momentum space:  see in file" << std::endl
                  << "    Helicity modules \\Rho_s:   " << dynamicMeasure.obs_mean_rho_s
                  << "    err: " << dynamicMeasure.obs_err_rho_s << std::endl;
        std::cout.precision(-1);
    }

    std::cout << std::endl;
    std::cout << "  Time Cost:      " << minute << " min " << sec << " s" << std::endl;

    std::cout << "==============================================================================" << std::endl;
}

void detQMC::output_Stats_eqtime(const std::string &filename, bool bool_Append) {
    if (bool_measure_eqtime) {
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
                << std::setw(15) << eqtimeMeasure.obs_mean_eqtime["DoubleOccu"]
                << std::setw(15) << eqtimeMeasure.obs_mean_eqtime["KineticEnergy"]
                << std::setw(15) << eqtimeMeasure.obs_mean_eqtime["StructFactor"]
                << std::setw(15) << eqtimeMeasure.obs_mean_eqtime["MomentumDist"]
                << std::setw(15) << eqtimeMeasure.obs_mean_eqtime["localSpinCorr"]
                << std::setw(15) << eqtimeMeasure.obs_err_eqtime["DoubleOccu"]
                << std::setw(15) << eqtimeMeasure.obs_err_eqtime["KineticEnergy"]
                << std::setw(15) << eqtimeMeasure.obs_err_eqtime["StructFactor"]
                << std::setw(15) << eqtimeMeasure.obs_err_eqtime["MomentumDist"]
                << std::setw(15) << eqtimeMeasure.obs_err_eqtime["localSpinCorr"]
                << std::setw(15) << eqtimeMeasure.q(0)
                << std::setw(15) << eqtimeMeasure.q(1)
                << std::endl;
        outfile.close();
        std::cout << "  Equal-time data has been written into file: " << filename << std::endl;
        if (! bool_measure_dynamic) {
            std::cout << "==============================================================================" << std::endl << std::endl;
        }
    }
}

void detQMC::output_Stats_dynamic(const std::string& filename, bool bool_Append) {
    if (bool_measure_dynamic) {
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
                    << std::setw(15) << dynamicMeasure.obs_mean_g_kt[l-1]
                    << std::setw(15) << dynamicMeasure.obs_err_g_kt[l-1]
                    << std::setw(15) << dynamicMeasure.obs_err_g_kt[l-1] / dynamicMeasure.obs_mean_g_kt[l-1]
                    << std::endl;
        }

        outfile << std::setw(15) << dynamicMeasure.obs_mean_rho_s
                << std::setw(15) << dynamicMeasure.obs_err_rho_s
                << std::setw(15) << dynamicMeasure.obs_err_rho_s / dynamicMeasure.obs_mean_rho_s
                << std::endl;

        outfile.close();
        std::cout << "  Dynamic data has been written into file: " << filename << std::endl;
        std::cout << "==============================================================================" << std::endl << std::endl;
    }
}

detQMC::~detQMC() {
    std::cout << std::endl << "The simulation was done :)" << std::endl;
}
