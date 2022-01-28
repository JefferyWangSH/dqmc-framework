#include "detqmc.h"
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

void Simulation::DetQMC::set_controlling_params(bool is_warm_up, bool is_eqtime_measure, bool is_dynamic_measure, bool is_checkerboard) {
    this->is_warm_up = is_warm_up;
    this->is_eqtime_measure = is_eqtime_measure;
    this->is_dynamic_measure = is_dynamic_measure;
    if (!this->hubbard) { this->hubbard = std::unique_ptr<Model::Hubbard>(new Model::Hubbard()); }
    this->hubbard->set_bool_params(is_eqtime_measure, is_dynamic_measure, is_checkerboard);
}

void Simulation::DetQMC::set_observable_list(const std::vector<std::string> &obs_list) {
    if (this->obs_list) { this->obs_list.reset(); }
    this->obs_list = std::unique_ptr<std::vector<std::string>>(new std::vector<std::string>(obs_list));
}

void Simulation::DetQMC::set_aux_field_configs(const std::string &config_file) {
    if (this->config_file) { this->config_file.reset(); }
    this->config_file = std::unique_ptr<std::string>(new std::string(config_file));
}

void Simulation::DetQMC::set_lattice_momentum(double qx, double qy) {
    this->q = (Eigen::VectorXd(2) << qx, qy).finished();
    // useful when scanning the momentum space
    if (this->measure) { this->measure->set_lattice_momentum(M_PI*this->q); }
}

void Simulation::DetQMC::read_configs_from_file(const std::string &config_file) {
    // caution: model params should be set up ahead
    std::ifstream infile;
    infile.open(config_file, std::ios::in);
    if (!infile.is_open()) {
        std::cerr << " Fail to read configs from file " + config_file + " ." << std::endl;
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
    if ( (lt != this->hubbard->lt) || (ls != this->hubbard->ls) ) {
        std::cerr << " Inconsistency between model settings and input configs (lt or ll). " << std::endl;
        exit(1);
    }

    // read in configs of aux field
    int time_point, space_point;
    while(getline(infile, line)) {
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
        time_point = boost::lexical_cast<int>(data[0]);
        space_point = boost::lexical_cast<int>(data[1]);
        (*(this->hubbard->s))(space_point, time_point) = boost::lexical_cast<double>(data[2]);
    }
    infile.close();
}

void Simulation::DetQMC::initial_hubbard_with_input_configs() {
    // allocate memory
    this->hubbard->allocate();

    // initialize checkerboard
    this->hubbard->checkerboard->init_from_model(*this->hubbard);
    this->hubbard->is_checkerboard = this->hubbard->checkerboard->is_checker_board();

    // initialize configs of aux fields from input file 
    this->read_configs_from_file(*this->config_file);

    // refresh greens and svd stacks for input configs
    this->hubbard->init_stacks();
    this->hubbard->config_sign = (this->hubbard->green_tt_up->determinant() * this->hubbard->green_tt_dn->determinant() >= 0)? +1.0 : -1.0;
}

void Simulation::DetQMC::initial() {

    // initialize model, with configs of aux fields set to random
    if (this->is_warm_up) { this->hubbard->initial(); }
    // initialize model with input configs
    else if (this->config_file) { this->initial_hubbard_with_input_configs(); }
    else { 
        this->hubbard->initial();
        std::cerr << " Caution: no input configurations and not warm up, the simulation results may be incorrect. " << std::endl; 
    }

    // initialize measure class
    if (this->measure) { this->measure.reset(); }
    if (this->obs_list) {
        this->measure = std::unique_ptr<Measure::Measure>(new Measure::Measure());
        this->measure->set_size_of_bin(this->nbin);
        this->measure->set_observable_list(*this->obs_list);
        this->measure->set_lattice_momentum(M_PI*this->q);
        this->measure->initial(*this->hubbard);
    }

    // consistency check
    if (this->measure) {
        if ((this->measure->_container.is_eqtime_measure() != this->is_eqtime_measure) || 
            (this->measure->_container.is_dynamic_measure() != this->is_dynamic_measure)) {
            std::cerr << " Consistency check failed, check the input params of observable measurments. " << std::endl;
            exit(1);
        }
    }
    else if (this->is_eqtime_measure || this->is_dynamic_measure) {
        std::cerr << " Consistency check failed, check the input params of observable measurments. " << std::endl;
        exit(1);
    }
    if (!this->is_eqtime_measure && !this->is_dynamic_measure) { this->measure.reset(); }

}

void Simulation::DetQMC::sweep_forth_and_back(bool is_eqtime_measure, bool is_dynamic_measure) const {
    // sweep forth from 0 to beta
    if (!is_dynamic_measure) {
        this->hubbard->sweep_0_to_beta();
    }
    else {
        this->hubbard->sweep_0_to_beta_dynamic();
        this->measure->dynamic_measure(*this->hubbard);
    }
    if (is_eqtime_measure) {
        this->measure->eqtime_measure(*this->hubbard);
    }

    // sweep back from beta to 0
    this->hubbard->sweep_beta_to_0();
    // TODO: hubb.sweep_beta_to_0_displaced
    if (is_eqtime_measure) {
        this->measure->eqtime_measure(*this->hubbard);
    }
}

void Simulation::DetQMC::run(bool show_running_process) {
    assert( this->hubbard );

    // clear previous temporary data in case
    if (this->measure) { this->measure->clear_temporary(); }

    // record current time
    this->begin_t = std::chrono::steady_clock::now();

    // thermalization process
    if (this->is_warm_up) {
        // progress bar
        progresscpp::ProgressBar progress_bar(this->nwarm/2, 40, '#', '-');

        for (int nwm = 1; nwm <= this->nwarm/2; ++nwm) {
            this->sweep_forth_and_back(false, false);
            ++progress_bar;

            if ( nwm % 10 == 0 && show_running_process ) {
                std::cout << " Warm-up progress:   "; 
                progress_bar.display(); 
            }
        }
        if (show_running_process) {
            std::cout << " Warm-up progress:   ";
            progress_bar.done();
        }
    }

    if (this->measure) {
        // measuring process
        progresscpp::ProgressBar progress_bar(this->nbin * this->nsweep / 2, 40, '#', '-');

        for (int bin = 0; bin < this->nbin; ++bin) {
            for (int nsw = 1; nsw <= this->nsweep/2; ++nsw) {
                this->sweep_forth_and_back(this->is_eqtime_measure, this->is_dynamic_measure);
                ++progress_bar;

                if ( nsw % 10 == 0 && show_running_process ) {
                    std::cout << " Measuring progress: ";
                    progress_bar.display();
                }
            }

            // analyse collected data
            if (this->measure) {
                this->measure->normalize_stats();
                this->measure->write_stats_to_bins(bin);
                this->measure->clear_temporary();
            }

            // avoid correlation between bins
            for (int nbtw = 0; nbtw < this->n_between_bins; ++nbtw) {
                this->sweep_forth_and_back(false, false);
            }
        }
        if (show_running_process) {
            std::cout << " Measuring progress: ";
            progress_bar.done();
        }
    }
    // record ending time
    this->end_t = std::chrono::steady_clock::now();
}

void Simulation::DetQMC::analyse_stats() const {
    if (this->measure) {
        this->measure->analyse_stats();
    }
}

