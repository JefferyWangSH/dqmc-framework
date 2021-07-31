#include "detQMC.h"
#include "ProgressBar.hpp"
#include <cmath>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

void detQMC::set_Model_Params(int _ll, int _lt, double _beta, double _t, double _Uint, double _mu, int _nwrap) {
    Hubbard newHubb(_ll, _lt, _beta, _t, _Uint, _mu, nwrap);
    hubb = newHubb;
    this->nwrap = _nwrap;
}

void detQMC::set_MC_Params(int _nwarm, int _nbin, int _nsweep, int _nBetweenBins) {
    this->nwarm = _nwarm;
    this->nsweep = _nsweep;
    this->nBetweenBins = _nBetweenBins;
    this->nbin = _nbin;

    eqtimeMeasure.resize(_nbin);
    dynamicMeasure.resize(_nbin);
}

void detQMC::set_bool_Params(bool _bool_warm_up, bool _bool_measure_eqtime, bool _bool_measure_dynamic) {
    this->bool_warm_up = _bool_warm_up;
    this->bool_measure_eqtime = _bool_measure_eqtime;
    this->bool_measure_dynamic = _bool_measure_dynamic;
}

void detQMC::set_Momentum_q(double qx, double qy) {
    q = (Eigen::VectorXd(2) << qx, qy).finished();
    eqtimeMeasure.q = (Eigen::VectorXd(2) << qx, qy).finished();
    dynamicMeasure.q = (Eigen::VectorXd(2) << qx, qy).finished();
}

void detQMC::read_Aux_Field_Configs(const std::string &filename) {
    // model params should be set up ahead

    std::ifstream infile;
    infile.open(filename, std::ios::in);

    if (!infile.is_open()) {
        std::cerr << "fail to open file " + filename + " !" << std::endl;
        exit(1);
    }

    std::string line;
    int lt = 0;
    int ls = 0;
    while(getline(infile, line)) {
        std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));

        int l = boost::lexical_cast<int>(data[0]);
        int i = boost::lexical_cast<int>(data[1]);
        hubb.s(i, l) = boost::lexical_cast<double>(data[2]);
        lt = std::max(lt, l);
        ls = std::max(ls, i);
    }
    assert(lt + 1 == hubb.lt);
    assert(ls + 1 == hubb.ls);
    infile.close();

    /* initial greens and svd stacks for old configs */
    hubb.stackLeftU = new SvdStack(hubb.ls, hubb.lt);
    hubb.stackLeftD = new SvdStack(hubb.ls, hubb.lt);
    hubb.stackRightU = new SvdStack(hubb.ls, hubb.lt);
    hubb.stackRightD = new SvdStack(hubb.ls, hubb.lt);
    hubb.initStacks(nwrap);
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

        for (int bin = 0; bin < nbin; ++bin) {
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

    std::cout << std::endl;
    std::cout << "  Maximum of wrap error (equal-time):     " << hubb.max_wrap_error_equal << std::endl
              << "  Maximum of wrap error (time-displaced): " << hubb.max_wrap_error_displaced << std::endl;
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
    // TODO: hubb.sweep_beta_to_0_displaced
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
                  << "    Double Occupancy:        " << eqtimeMeasure.obs_mean_eqtime["double_occupancy"]
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["double_occupancy"] << std::endl
                  << "    Kinetic Energy:          " << eqtimeMeasure.obs_mean_eqtime["kinetic_energy"]
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["kinetic_energy"] << std::endl
                  << "    Momentum Distribution:   " << eqtimeMeasure.obs_mean_eqtime["momentum_distribution"]
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["momentum_distribution"] << std::endl
                  << "    Local Spin Correlation:  " << eqtimeMeasure.obs_mean_eqtime["local_spin_correlation"]
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["local_spin_correlation"] << std::endl
                  << "    Structure Factor:        " << eqtimeMeasure.obs_mean_eqtime["structure_factor"]
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["structure_factor"] << std::endl
                  << "    Average Sign (abs):      " << abs(eqtimeMeasure.obs_mean_eqtime["average_sign"])
                  << "    err: " << eqtimeMeasure.obs_err_eqtime["average_sign"] << std::endl;
        std::cout.precision(-1);
    }

    if (bool_measure_dynamic) {
        std::cout.precision(8);
        std::cout << std::endl;
        std::cout << "  Time-displaced Measurements: " << std::endl
                  << "    Dynamical correlation in momentum space:  see in file" << std::endl
                  << "    Correlation G(k, beta/2):   " << dynamicMeasure.obs_mean_g_kt[ceil(hubb.lt/2.0)]
                  << "    err: " << dynamicMeasure.obs_err_g_kt[ceil(hubb.lt/2.0)] << std::endl
                  << "    Helicity modules \\Rho_s:   " << dynamicMeasure.obs_mean_rho_s
                  << "    err: " << dynamicMeasure.obs_err_rho_s << std::endl
                  << "    Average Sign (abs):         " << abs(dynamicMeasure.obs_mean_sign)
                  << "    err: " << dynamicMeasure.obs_err_sign << std::endl;
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
                << std::setw(15) << eqtimeMeasure.obs_mean_eqtime["double_occupancy"]
                << std::setw(15) << eqtimeMeasure.obs_mean_eqtime["kinetic_energy"]
                << std::setw(15) << eqtimeMeasure.obs_mean_eqtime["structure_factor"]
                << std::setw(15) << eqtimeMeasure.obs_mean_eqtime["momentum_distribution"]
                << std::setw(15) << eqtimeMeasure.obs_mean_eqtime["local_spin_correlation"]
                << std::setw(15) << eqtimeMeasure.obs_err_eqtime["double_occupancy"]
                << std::setw(15) << eqtimeMeasure.obs_err_eqtime["kinetic_energy"]
                << std::setw(15) << eqtimeMeasure.obs_err_eqtime["structure_factor"]
                << std::setw(15) << eqtimeMeasure.obs_err_eqtime["momentum_distribution"]
                << std::setw(15) << eqtimeMeasure.obs_err_eqtime["local_spin_correlation"]
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

void detQMC::output_Aux_Field_Configs(const std::string &filename) {
    std::ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::trunc);

    outfile << std::setiosflags(std::ios::right);
    for (int l = 0; l < hubb.lt; ++l) {
        for (int i = 0; i < hubb.ls; ++i) {
            outfile << std::setw(15) << l
                    << std::setw(15) << i
                    << std::setw(15) << hubb.s(i, l)
                    << std::endl;
        }
    }
}

detQMC::~detQMC() {
    std::cout << std::endl << "The simulation was done :)" << std::endl;
}
