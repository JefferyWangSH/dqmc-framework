#include "DetQMC.h"
#include "Hubbard.h"
#include "SvdStack.h"
#include "EqtimeMeasure.h"
#include "DynamicMeasure.h"
#include "ProgressBar.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


void Simulation::DetQMC::set_model_params(int _ll, int _lt, double _beta, double _t, double _uint, double _mu, int _nwrap, bool _is_checkerboard) {
    this->nwrap = _nwrap;
    if (this->hubb) {
        delete this->hubb;
        this->hubb = new Model::Hubbard(_ll, _lt, _beta, _t, _uint, _mu, _nwrap, _is_checkerboard);
    }
    else {
        this->hubb = new Model::Hubbard(_ll, _lt, _beta, _t, _uint, _mu, _nwrap, _is_checkerboard);
    }
}

void Simulation::DetQMC::set_Monte_Carlo_params(int _nwarm, int _nbin, int _nsweep, int _n_between_bins) {
    this->nwarm = _nwarm;
    this->nsweep = _nsweep;
    this->n_between_bins = _n_between_bins;
    this->nbin = _nbin;
}

void Simulation::DetQMC::set_controlling_params(bool _bool_warm_up, bool _bool_measure_eqtime, bool _bool_measure_dynamic) {
    this->bool_warm_up = _bool_warm_up;
    this->bool_measure_eqtime = _bool_measure_eqtime;
    this->bool_measure_dynamic = _bool_measure_dynamic;
}

void Simulation::DetQMC::set_lattice_momentum(double qx, double qy) {
    this->q = (Eigen::VectorXd(2) << qx, qy).finished();
    if (this->EqtimeMeasure) { this->EqtimeMeasure->q = M_PI * this->q; }
    if (this->DynamicMeasure) { this->DynamicMeasure->q = M_PI * this->q; }
}

void Simulation::DetQMC::read_aux_field_configs(const std::string &filename) {
    // Caution: Model params should be set up ahead
    std::ifstream infile;
    infile.open(filename, std::ios::in);

    if (!infile.is_open()) {
        std::cerr << "fail to open file " + filename + " !" << std::endl;
        exit(1);
    }

    int lt = 0, ls = 0;
    std::string line;
    while(getline(infile, line)) {
        std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));

        int l = boost::lexical_cast<int>(data[0]);
        int i = boost::lexical_cast<int>(data[1]);
        this->hubb->s(i, l) = boost::lexical_cast<double>(data[2]);
        lt = std::max(lt, l);
        ls = std::max(ls, i);
    }
    assert( lt + 1 == this->hubb->lt );
    assert( ls + 1 == this->hubb->ls );
    infile.close();

    /* initial greens and svd stacks for input configs */
    this->hubb->init_stacks(this->nwrap);
}

void Simulation::DetQMC::init_measure() {
    // initialize bins of observables
    // for equal-time measurements
    if (this->bool_measure_eqtime && this->EqtimeMeasure) {
        delete this->EqtimeMeasure;
        this->EqtimeMeasure = new Measure::EqtimeMeasure(this->nbin);
        this->EqtimeMeasure->initial(*this->hubb);
        this->EqtimeMeasure->q = M_PI * this->q;
    }
    else if (this->bool_measure_eqtime && !this->EqtimeMeasure) {
        this->EqtimeMeasure = new Measure::EqtimeMeasure(this->nbin);
        this->EqtimeMeasure->initial(*this->hubb);
        this->EqtimeMeasure->q = M_PI * this->q;
    }
    else { this->EqtimeMeasure = nullptr; }
    
    // for dynamical measurements
    if (this->bool_measure_dynamic && this->DynamicMeasure) {
        delete this->DynamicMeasure;
        this->DynamicMeasure = new Measure::DynamicMeasure(this->nbin);
        this->DynamicMeasure->initial(*this->hubb);
        this->DynamicMeasure->q = M_PI * this->q;
    }
    else if (this->bool_measure_dynamic && !this->DynamicMeasure) {
        this->DynamicMeasure = new Measure::DynamicMeasure(this->nbin);
        this->DynamicMeasure->initial(*this->hubb);
        this->DynamicMeasure->q = M_PI * this->q;
    }
    else { this->DynamicMeasure = nullptr; }
}

void Simulation::DetQMC::run_QMC(bool bool_display_process) {
    assert( this->hubb );

    // clear data
    if (this->bool_measure_eqtime && this->EqtimeMeasure) {
        this->EqtimeMeasure->clear_temporary();
    }
    if (this->bool_measure_dynamic && this->DynamicMeasure) {
        this->DynamicMeasure->clear_temporary(*hubb);
    }

    // record current time
    this->begin_t = std::chrono::steady_clock::now();

    // thermalization process
    if (this->bool_warm_up) {
        // progress bar
        progresscpp::ProgressBar progressBar(nwarm/2, 40, '#', '-');

        for (int nwm = 1; nwm <= nwarm/2; ++nwm) {
            this->sweep_back_and_forth(false, false);
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

    if (this->bool_measure_eqtime || this->bool_measure_dynamic) {
        // measuring process
        progresscpp::ProgressBar progressBar(nbin * nsweep / 2, 40, '#', '-');

        for (int bin = 0; bin < nbin; ++bin) {
            for (int nsw = 1; nsw <= nsweep/2; ++nsw) {
                this->sweep_back_and_forth(this->bool_measure_eqtime, this->bool_measure_dynamic);
                ++progressBar;

                if ( nsw % 10 == 0 && bool_display_process) {
                    std::cout << "Measuring progress: ";
                    progressBar.display();
                }
            }

            // analyse statistical data
            if (this->bool_measure_eqtime && this->EqtimeMeasure) {
                this->EqtimeMeasure->normalize_stats(*this->hubb);
                this->EqtimeMeasure->write_stats_to_bins(bin);
                this->EqtimeMeasure->clear_temporary();
            }

            if (this->bool_measure_dynamic && this->DynamicMeasure) {
                this->DynamicMeasure->normalize_stats(*this->hubb);
                this->DynamicMeasure->write_stats_to_bins(bin, *this->hubb);
                this->DynamicMeasure->clear_temporary(*this->hubb);
            }

            // avoid correlation between bins
            for (int n_bw = 0; n_bw < n_between_bins; ++n_bw) {
                this->sweep_back_and_forth(false, false);
            }
        }

        if (bool_display_process) {
            std::cout << "Measuring progress: ";
            progressBar.done();
        }
    }

    std::cout << std::endl;
    std::cout << "  Maximum of wrap error (equal-time):     " << this->hubb->max_wrap_error_equal << std::endl
              << "  Maximum of wrap error (time-displaced): " << this->hubb->max_wrap_error_displaced << std::endl;
    this->end_t = std::chrono::steady_clock::now();
}

void Simulation::DetQMC::sweep_back_and_forth(bool bool_eqtime, bool bool_dynamic) const {

    // sweep forth from 0 to beta
    if (!bool_dynamic) {
        this->hubb->sweep_0_to_beta(this->nwrap);
    }
    else {
        this->hubb->sweep_0_to_beta_displaced(this->nwrap);
        this->DynamicMeasure->time_displaced_measure(*this->hubb);
    }
    if (bool_eqtime) {
        this->EqtimeMeasure->equal_time_measure(*this->hubb);
    }

    // sweep back from beta to 0
    this->hubb->sweep_beta_to_0(this->nwrap);
    // TODO: hubb.sweep_beta_to_0_displaced
    if (bool_eqtime) {
        this->EqtimeMeasure->equal_time_measure(*this->hubb);
    }
}

void Simulation::DetQMC::analyse_stats() const {
    if (this->bool_measure_eqtime) {
        this->EqtimeMeasure->analyse_stats(*this->hubb);
    }
    if (this->bool_measure_dynamic) {
        this->DynamicMeasure->analyse_stats(*this->hubb);
    }
}

void Simulation::DetQMC::print_params() const{
    std::cout << std::endl;
    std::cout << "==============================================================================" << std::endl;
    std::cout << "  Simulation Parameters: " << std::endl
    << "    ll:     " << this->hubb->ll << std::endl
    << "    lt:     " << this->hubb->lt << std::endl
    << "    beta:   " << this->hubb->beta << std::endl
    << "    U/t:    " << this->hubb->Uint / this->hubb->t << std::endl
    << "    mu:     " << this->hubb->mu << std::endl
    << "    q:      " << this->q(0) << " pi, "<< this->q(1) << " pi" << std::endl
    << "    nwrap:  " << this->nwrap << std::endl;
    std::cout << "==============================================================================" << std::endl;
}

void Simulation::DetQMC::print_stats() const {
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(this->end_t - this->begin_t).count();
    const int minute = std::floor((double)time / 1000 / 60);
    const double sec = (double)time / 1000 - 60 * minute;

    if (this->bool_measure_eqtime) {
        std::cout.precision(8);
        std::cout << std::endl;
        std::cout << "  Equal-time Measurements: " << std::endl
                  << "    Double Occupancy:        " << this->EqtimeMeasure->double_occu.mean_value()
                  << "    err: " << this->EqtimeMeasure->double_occu.error_bar() << std::endl
                  << "    Kinetic Energy:          " << this->EqtimeMeasure->kinetic_energy.mean_value()
                  << "    err: " << this->EqtimeMeasure->kinetic_energy.error_bar() << std::endl
                  << "    Electron Density:        " << this->EqtimeMeasure->electron_density.mean_value()
                  << "    err: " << this->EqtimeMeasure->electron_density.error_bar() << std::endl
                  << "    Local Spin Correlation:  " << this->EqtimeMeasure->local_corr.mean_value()
                  << "    err: " << this->EqtimeMeasure->local_corr.error_bar() << std::endl
                  << "    AFM Structure Factor:    " << this->EqtimeMeasure->AFM_factor.mean_value()
                  << "    err: " << this->EqtimeMeasure->AFM_factor.error_bar() << std::endl
                  << "    Average Sign (abs):      " << abs(this->EqtimeMeasure->sign.mean_value())
                  << "    err: " << this->EqtimeMeasure->sign.error_bar() << std::endl;
        std::cout.precision(-1);
    }

    if (this->bool_measure_dynamic) {
        std::cout.precision(8);
        std::cout << std::endl;
        std::cout << "  Time-displaced Measurements: " << std::endl
                  << "    Dynamical green's functions of imaginary time:  see in file" << std::endl
                  << "    Green's function G(k, beta/2):   " << this->DynamicMeasure->matsubara_greens[ceil(this->hubb->lt/2.0)].mean_value()
                  << "    err: " << this->DynamicMeasure->matsubara_greens[ceil(this->hubb->lt/2.0)].error_bar() << std::endl
                  << "    Superfluid density rho_s:   " << this->DynamicMeasure->superfluid_density.mean_value()
                  << "    err: " << this->DynamicMeasure->superfluid_density.error_bar() << std::endl
                  << "    Average Sign (abs):         " << abs(this->DynamicMeasure->sign.mean_value())
                  << "    err: " << this->DynamicMeasure->sign.error_bar() << std::endl;
        std::cout.precision(-1);
    }

    std::cout << std::endl;
    std::cout << "  Time Cost:      " << minute << " min " << sec << " s" << std::endl;

    std::cout << "==============================================================================" << std::endl;
}

void Simulation::DetQMC::file_output_tau(const std::string &filename) const {
    std::ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::trunc);

    outfile << std::setiosflags(std::ios::right)
            << std::setw(7) << this->hubb->lt
            << std::setw(7) << this->hubb->beta << std::endl;
    for (int l = 0; l < this->hubb->lt; ++l){
        outfile << std::setw(15) << l * this->hubb->dtau << std::endl;
    }
    outfile.close();
}

void Simulation::DetQMC::bin_output_greens(const std::string &filename) const{
    if (this->bool_measure_dynamic) {
        std::ofstream outfile;
        std::ifstream test_exist(filename);
        if (!test_exist.good()) {
            outfile.open(filename, std::ios::out | std::ios::app);
            outfile << std::setiosflags(std::ios::right) << std::setw(10) << this->nbin << std::endl;
            outfile.close();
        }
        else { test_exist.close(); }

        outfile.open(filename, std::ios::out | std::ios::app);
        outfile.precision(15);
        for (int bin = 0; bin < this->nbin; ++bin) {
            outfile << std::setw(20) << bin << std::endl;
            for (int l = 0; l < this->hubb->lt; ++l) {
                const int tau = (l - 1 + this->hubb->lt) % this->hubb->lt;
                outfile << std::setw(20) << this->DynamicMeasure->matsubara_greens[tau].bin_data()[bin] << std::endl;
            }
        }
        outfile.close();
    }
}

void Simulation::DetQMC::bin_output_LDOS(const std::string &filename) const {
    if (this->bool_measure_dynamic) {
        std::ofstream outfile;
        std::ifstream test_exist(filename);
        if (!test_exist.good()) {
            outfile.open(filename, std::ios::out | std::ios::app);
            outfile << std::setiosflags(std::ios::right) << std::setw(10) << this->nbin << std::endl;
            outfile.close();
        }
        else { test_exist.close(); }

        outfile.open(filename, std::ios::out | std::ios::app);
        outfile.precision(15);
        for (int bin = 0; bin < this->nbin; ++bin) {
            outfile << std::setw(20) << bin << std::endl;
            for (int l = 0; l < this->hubb->lt; ++l) {
                const int tau = (l - 1 + this->hubb->lt) % this->hubb->lt;
                // TODO: LDOS measurements (linear) can be accelerated!
                outfile << std::setw(20)
                        << this->DynamicMeasure->density_of_states[tau].bin_data()[bin]
                        << std::endl;
            }
        }
        outfile.close();
    }
}

void Simulation::DetQMC::file_output_eqtime_stats(const std::string &filename) const {
    if (this->bool_measure_eqtime) {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out | std::ios::trunc);

        outfile << std::setiosflags(std::ios::right)
                << std::setw(15) << this->hubb->Uint / this->hubb->t
                << std::setw(15) << this->hubb->beta
                << std::setw(15) << this->EqtimeMeasure->double_occu.mean_value()
                << std::setw(15) << this->EqtimeMeasure->double_occu.error_bar()
                << std::setw(15) << this->EqtimeMeasure->kinetic_energy.mean_value()
                << std::setw(15) << this->EqtimeMeasure->kinetic_energy.error_bar()
                << std::setw(15) << this->EqtimeMeasure->AFM_factor.mean_value()
                << std::setw(15) << this->EqtimeMeasure->AFM_factor.error_bar()
                << std::setw(15) << this->EqtimeMeasure->electron_density.mean_value()
                << std::setw(15) << this->EqtimeMeasure->electron_density.error_bar()
                << std::setw(15) << this->EqtimeMeasure->local_corr.mean_value()
                << std::setw(15) << this->EqtimeMeasure->local_corr.error_bar()
                << std::setw(15) << this->EqtimeMeasure->q(0)
                << std::setw(15) << this->EqtimeMeasure->q(1)
                << std::endl;
        outfile.close();
        std::cout << "  Equal-time data has been written into file: " << filename << std::endl;
        if (! this->bool_measure_dynamic) {
            std::cout << "==============================================================================" << std::endl << std::endl;
        }
    }
}

void Simulation::DetQMC::file_output_dynamic_stats(const std::string& filename) const{
    if (this->bool_measure_dynamic) {
        std::ofstream outfile;
        outfile.open(filename, std::ios::out | std::ios::trunc);

        outfile << std::setiosflags(std::ios::right)
                << "Momentum k: " << this->q(0) << " pi, "<< this->q(1) << " pi" << std::endl;

        for (int l = 0; l < this->hubb->lt; ++l) {
            const int tau = (l - 1 + this->hubb->lt) % this->hubb->lt;
            outfile << std::setw(15) << l
                    << std::setw(15) << this->DynamicMeasure->matsubara_greens[tau].mean_value()
                    << std::setw(15) << this->DynamicMeasure->matsubara_greens[tau].error_bar()
                    << std::setw(15) << this->DynamicMeasure->matsubara_greens[tau].error_bar() / this->DynamicMeasure->matsubara_greens[tau].mean_value()
                    << std::endl;
        }

        outfile << std::setw(15) << this->DynamicMeasure->superfluid_density.mean_value()
                << std::setw(15) << this->DynamicMeasure->superfluid_density.error_bar()
                << std::setw(15) << this->DynamicMeasure->superfluid_density.error_bar() / this->DynamicMeasure->superfluid_density.mean_value()
                << std::endl;

        outfile.close();
        std::cout << "  Dynamic data has been written into file: " << filename << std::endl;
        std::cout << "==============================================================================" << std::endl << std::endl;
    }
}

void Simulation::DetQMC::file_output_cooper_corr(const std::string &filename) const {
    std::ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::trunc);

    for (int i = 0; i < this->EqtimeMeasure->cooper_corr.size(); ++i) {
        outfile << std::setw(15) << i
                << std::setw(15) << this->EqtimeMeasure->cooper_corr[i].mean_value()
                << std::setw(15) << this->EqtimeMeasure->cooper_corr[i].error_bar()
                << std::endl;
    }
    outfile.close();
}

void Simulation::DetQMC::file_output_aux_field_configs(const std::string &filename) const{
    std::ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::trunc);

    outfile << std::setiosflags(std::ios::right);
    for (int l = 0; l < this->hubb->lt; ++l) {
        for (int i = 0; i < this->hubb->ls; ++i) {
            outfile << std::setw(15) << l
                    << std::setw(15) << i
                    << std::setw(15) << this->hubb->s(i, l)
                    << std::endl;
        }
    }
    outfile.close();
}

Simulation::DetQMC::~DetQMC() {
    if (this->hubb) {
        delete this->hubb;
        this->hubb = nullptr;
    }
    if (this->EqtimeMeasure) {
        delete this->EqtimeMeasure;
        this->EqtimeMeasure = nullptr;
    }
    if (this->DynamicMeasure) {
        delete DynamicMeasure;
        this->DynamicMeasure = nullptr;
    }
    std::cout << std::endl << "The simulation was done :)" << std::endl;
}
