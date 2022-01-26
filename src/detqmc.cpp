#include "detqmc.h"
#include "hubbard.h"
#include "svd_stack.h"
#include "eqtime_measure.h"
#include "dynamic_measure.h"
#include "progress_bar.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


void Simulation::DetQMC::set_model_params(int ll, int lt, double beta, double t, double u_int, double mu, int nwrap) {
    if (!this->hubbard) { this->hubbard = std::unique_ptr<Model::Hubbard>(new Model::Hubbard()); }
    this->hubbard->set_model_params(ll, lt, beta, t, u_int, mu);
    this->hubbard->set_stabilization_pace(nwrap);
    this->nwrap = nwrap;
}

void Simulation::DetQMC::set_Monte_Carlo_params(int nwarm, int nbin, int nsweep, int n_between_bins) {
    this->nwarm = nwarm;
    this->nbin = nbin;
    this->nsweep = nsweep;
    this->n_between_bins = n_between_bins;
}

void Simulation::DetQMC::set_controlling_params(bool _bool_warm_up, bool _bool_measure_eqtime, bool _bool_measure_dynamic, bool is_checkerboard) {
    this->bool_warm_up = _bool_warm_up;
    this->bool_measure_eqtime = _bool_measure_eqtime;
    this->bool_measure_dynamic = _bool_measure_dynamic;
    if (!this->hubbard) { this->hubbard = std::unique_ptr<Model::Hubbard>(new Model::Hubbard()); }
    this->hubbard->set_bool_params(_bool_measure_eqtime, _bool_measure_dynamic, is_checkerboard);
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
        std::cerr << " Fail to read configs from file " + filename + " !" << std::endl;
        exit(1);
    }

    // temporary parmas
    std::string line;
    std::vector<std::string> data;
    
    // check the model parameters
    getline(infile, line);
    boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
    data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
    const int lt = boost::lexical_cast<int>(data[0]);
    const int ls = boost::lexical_cast<int>(data[1]); 
    assert( lt == this->hubbard->lt );
    assert( ls == this->hubbard->ls );

    // read in configs of aux field
    int time_point, space_point;
    while(getline(infile, line)) {
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
        time_point = boost::lexical_cast<int>(data[0]);
        space_point = boost::lexical_cast<int>(data[1]);
        (*this->hubbard->s)(space_point, time_point) = boost::lexical_cast<double>(data[2]);
    }
    infile.close();

    /* initial greens and svd stacks for input configs */
    this->hubbard->init_stacks();
    this->hubbard->config_sign = ((*this->hubbard->green_tt_up).determinant() * (*this->hubbard->green_tt_dn).determinant() >= 0)? +1.0 : -1.0;
}

void Simulation::DetQMC::initial() {
    // initialize model
    this->hubbard->initial();

    // initialize bins of observables
    // for equal-time measurements
    if (this->bool_measure_eqtime && this->EqtimeMeasure) {
        this->EqtimeMeasure.reset();
        this->EqtimeMeasure = std::unique_ptr<Measure::EqtimeMeasure>(new Measure::EqtimeMeasure(this->nbin));
        this->EqtimeMeasure->initial(*this->hubbard);
        this->EqtimeMeasure->q = M_PI * this->q;
    }
    else if (this->bool_measure_eqtime && !this->EqtimeMeasure) {
        this->EqtimeMeasure = std::unique_ptr<Measure::EqtimeMeasure>(new Measure::EqtimeMeasure(this->nbin));
        this->EqtimeMeasure->initial(*this->hubbard);
        this->EqtimeMeasure->q = M_PI * this->q;
    }
    else { this->EqtimeMeasure.reset(); }
    
    // for dynamical measurements
    if (this->bool_measure_dynamic && this->DynamicMeasure) {
        this->DynamicMeasure.reset();
        this->DynamicMeasure = std::unique_ptr<Measure::DynamicMeasure>(new Measure::DynamicMeasure(this->nbin));
        this->DynamicMeasure->initial(*this->hubbard);
        this->DynamicMeasure->q = M_PI * this->q;
    }
    else if (this->bool_measure_dynamic && !this->DynamicMeasure) {
        this->DynamicMeasure = std::unique_ptr<Measure::DynamicMeasure>(new Measure::DynamicMeasure(this->nbin));
        this->DynamicMeasure->initial(*this->hubbard);
        this->DynamicMeasure->q = M_PI * this->q;
    }
    else { this->DynamicMeasure.reset(); }
}

void Simulation::DetQMC::run(bool bool_display_process) {
    assert( this->hubbard );

    // clear data
    if (this->bool_measure_eqtime && this->EqtimeMeasure) {
        this->EqtimeMeasure->clear_temporary();
    }
    if (this->bool_measure_dynamic && this->DynamicMeasure) {
        this->DynamicMeasure->clear_temporary(*this->hubbard);
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
                std::cout << " Warm-up progress:   ";
                progressBar.display();
            }
        }

        if (bool_display_process) {
            std::cout << " Warm-up progress:   ";
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
                    std::cout << " Measuring progress: ";
                    progressBar.display();
                }
            }

            // analyse statistical data
            if (this->bool_measure_eqtime && this->EqtimeMeasure) {
                this->EqtimeMeasure->normalize_stats(*this->hubbard);
                this->EqtimeMeasure->write_stats_to_bins(bin);
                this->EqtimeMeasure->clear_temporary();
            }

            if (this->bool_measure_dynamic && this->DynamicMeasure) {
                this->DynamicMeasure->normalize_stats(*this->hubbard);
                this->DynamicMeasure->write_stats_to_bins(bin, *this->hubbard);
                this->DynamicMeasure->clear_temporary(*this->hubbard);
            }

            // avoid correlation between bins
            for (int n_bw = 0; n_bw < n_between_bins; ++n_bw) {
                this->sweep_back_and_forth(false, false);
            }
        }

        if (bool_display_process) {
            std::cout << " Measuring progress: ";
            progressBar.done();
        }
    }
    this->end_t = std::chrono::steady_clock::now();
}

void Simulation::DetQMC::sweep_back_and_forth(bool bool_eqtime, bool bool_dynamic) const {

    // sweep forth from 0 to beta
    if (!bool_dynamic) {
        this->hubbard->sweep_0_to_beta();
    }
    else {
        this->hubbard->sweep_0_to_beta_dynamic();
        this->DynamicMeasure->time_displaced_measure(*this->hubbard);
    }
    if (bool_eqtime) {
        this->EqtimeMeasure->equal_time_measure(*this->hubbard);
    }

    // sweep back from beta to 0
    this->hubbard->sweep_beta_to_0();
    // TODO: hubb.sweep_beta_to_0_displaced
    if (bool_eqtime) {
        this->EqtimeMeasure->equal_time_measure(*this->hubbard);
    }
}

void Simulation::DetQMC::analyse_stats() const {
    if (this->bool_measure_eqtime) {
        this->EqtimeMeasure->analyse_stats(*this->hubbard);
    }
    if (this->bool_measure_dynamic) {
        this->DynamicMeasure->analyse_stats(*this->hubbard);
    }
}
