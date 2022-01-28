#include "measure_container.h"
#include "measure_methods.h"
#include "hubbard.h"
#include "measure.h"

namespace Measure {

    std::vector<std::string> Container::list() const {
        if (this->_obs_list) { return *this->_obs_list; }
        else { return std::vector<std::string>(); }
    }

    std::vector<std::string> Container::eqtime_list() const {
        if (this->_obs_list_eqtime) { return *this->_obs_list_eqtime; }
        else { return std::vector<std::string>(); }
    }

    std::vector<std::string> Container::dynamic_list() const {
        if (this->_obs_list_dynamic) { return *this->_obs_list_dynamic; }
        else { return std::vector<std::string>(); }
    }

    std::vector<Observable<double>> Container::obs_list_double() const {
        if (this->_obs_double) { return *this->_obs_double; }
        else { return std::vector<Observable<double>>(); }
    }

    std::vector<Observable<Eigen::VectorXd>> Container::obs_list_vector() const {
        if (this->_obs_vector) { return *this->_obs_vector; }
        else { return std::vector<Observable<Eigen::VectorXd>>(); }
    }

    std::vector<Observable<Eigen::MatrixXd>> Container::obs_list_matrix() const {
        if (this->_obs_matrix) { return *this->_obs_matrix; }
        else { return std::vector<Observable<Eigen::MatrixXd>>(); }
    }

    bool Container::is_eqtime_measure() const {
        if(this->_obs_list_eqtime) { return true; }
        else {return false; }
    }

    bool Container::is_dynamic_measure() const {
        if(this->_obs_list_dynamic) { return true; }
        else {return false; }
    }

    bool Container::is_eqtime_obs(std::string obs) const {
        return (std::find(this->_supported_obs_eqtime.begin(), this->_supported_obs_eqtime.end(), obs) != this->_supported_obs_eqtime.end());
    }

    bool Container::is_dynamic_obs(std::string obs) const {
        return (std::find(this->_supported_obs_dynamic.begin(), this->_supported_obs_dynamic.end(), obs) != this->_supported_obs_dynamic.end());
    }

    void Container::read_list(const std::vector<std::string> &obs_list) {
        // preprocessing
        // remove redundant input
        std::vector<std::string> tmp_list = obs_list;
        std::sort(tmp_list.begin(), tmp_list.end());
        tmp_list.erase(unique(tmp_list.begin(),tmp_list.end()), tmp_list.end());

        // catch unsupported input observables and report errors
        for (auto obs : tmp_list) {
            if ( !this->is_eqtime_obs(obs) && !this->is_dynamic_obs(obs) ) {
                std::cerr << " Unsupported observable type " << obs << ". " << std::endl;
                exit(1);
            }
        }
        if (this->_obs_list) { this->_obs_list.reset(); }
        this->_obs_list = std::unique_ptr<std::vector<std::string>>(new std::vector<std::string>(tmp_list));
    }

    void Container::filter() {
        // release memory if allocated before
        if (this->_obs_list_eqtime) { this->_obs_list_eqtime.reset(); }
        if (this->_obs_list_dynamic) { this->_obs_list_dynamic.reset(); }

        for (auto obs : this->list()) {
            // catch equal-time observables
            if ( this->is_eqtime_obs(obs) ) {
                if (!this->_obs_list_eqtime) {
                    this->_obs_list_eqtime = std::unique_ptr<std::vector<std::string>>(new std::vector<std::string>());
                }
                this->_obs_list_eqtime->push_back(obs);
            }
            // catch dynamical observables
            if ( this->is_dynamic_obs(obs) ) {
                if (!this->_obs_list_dynamic) {
                    this->_obs_list_dynamic = std::unique_ptr<std::vector<std::string>>(new std::vector<std::string>());
                }
                this->_obs_list_dynamic->push_back(obs);
            }
        }
    }

    void Container::allocate() {
        // release pre-allocated memory
        if (this->_obs_double) { this->_obs_double.reset(); }
        if (this->_obs_vector) { this->_obs_vector.reset(); }
        if (this->_obs_matrix) { this->_obs_matrix.reset(); }
        if (this->_sign_eqtime) { this->_sign_eqtime.reset(); }
        if (this->_sign_dynamic) { this->_sign_dynamic.reset(); }

        // allocate for measurements of equal-time observables
        for (auto obs : this->eqtime_list()) {
            if (obs == "filling_number") {
                if (!this->_obs_double) {
                    this->_obs_double = std::unique_ptr<std::vector<Observable<double>>>(new std::vector<Observable<double>>());
                }
                // set up name and method for specific physical observable
                // with bin size and zero element unassigned
                Observable<double> filling_number;
                filling_number.set_observable_name("filling_number");
                filling_number.add_method(Methods::measure_filling_number);
                this->_obs_double->push_back(filling_number);
            }

            if (obs == "double_occupancy") {
                if (!this->_obs_double) {
                    this->_obs_double = std::unique_ptr<std::vector<Observable<double>>>(new std::vector<Observable<double>>());
                }
                Observable<double> double_occupancy;
                double_occupancy.set_observable_name("double_occupancy");
                double_occupancy.add_method(Methods::measure_double_occupancy);
                this->_obs_double->push_back(double_occupancy);
            }

            if (obs == "kinetic_energy") {
                if (!this->_obs_double) {
                    this->_obs_double = std::unique_ptr<std::vector<Observable<double>>>(new std::vector<Observable<double>>());
                }
                Observable<double> kinetic_energy;
                kinetic_energy.set_observable_name("kinetic_energy");
                kinetic_energy.add_method(Methods::measure_kinetic_energy);
                this->_obs_double->push_back(kinetic_energy);
            }

            if (obs == "momentum_distribution") {
                if (!this->_obs_double) {
                    this->_obs_double = std::unique_ptr<std::vector<Observable<double>>>(new std::vector<Observable<double>>());
                }
                Observable<double> momentum_dist;
                momentum_dist.set_observable_name("momentum_distribution");
                momentum_dist.add_method(Methods::measure_momentum_distribution);
                this->_obs_double->push_back(momentum_dist);
            }

            if (obs == "local_spin_corr") {
                if (!this->_obs_double) {
                    this->_obs_double = std::unique_ptr<std::vector<Observable<double>>>(new std::vector<Observable<double>>());
                }
                Observable<double> local_spin_corr;
                local_spin_corr.set_observable_name("local_spin_corr");
                local_spin_corr.add_method(Methods::measure_local_spin_corr);
                this->_obs_double->push_back(local_spin_corr);
            }

            if (obs == "spin_density_structure_factor") {
                if (!this->_obs_double) {
                    this->_obs_double = std::unique_ptr<std::vector<Observable<double>>>(new std::vector<Observable<double>>());
                }
                Observable<double> sdw_factor;
                sdw_factor.set_observable_name("spin_density_structure_factor");
                sdw_factor.add_method(Methods::measure_spin_density_structure_factor);
                this->_obs_double->push_back(sdw_factor);
            }

            if (obs == "charge_density_structure_factor") {
                if (!this->_obs_double) {
                    this->_obs_double = std::unique_ptr<std::vector<Observable<double>>>(new std::vector<Observable<double>>());
                }
                Observable<double> cdw_factor;
                cdw_factor.set_observable_name("charge_density_structure_factor");
                cdw_factor.add_method(Methods::measure_charge_density_structure_factor);
                this->_obs_double->push_back(cdw_factor);
            }

            // adding new methods here
        }


        // allocate for measurements of dynamic observables
        for (auto obs : this->dynamic_list()) {
            if (obs == "matsubara_greens") {
                if (!this->_obs_vector) {
                    this->_obs_vector = std::unique_ptr<std::vector<Observable<Eigen::VectorXd>>>(new std::vector<Observable<Eigen::VectorXd>>());
                }
                Observable<Eigen::VectorXd> matsubara_greens;
                matsubara_greens.set_observable_name("matsubara_greens");
                matsubara_greens.add_method(Methods::measure_matsubara_greens);
                this->_obs_vector->push_back(matsubara_greens);
            }

            if (obs == "density_of_states") {
                if (!this->_obs_vector) {
                    this->_obs_vector = std::unique_ptr<std::vector<Observable<Eigen::VectorXd>>>(new std::vector<Observable<Eigen::VectorXd>>());
                }
                Observable<Eigen::VectorXd> density_of_states;
                density_of_states.set_observable_name("density_of_states");
                density_of_states.add_method(Methods::measure_density_of_states);
                this->_obs_vector->push_back(density_of_states);
            }

            if (obs == "superfluid_stiffness") {
                if (!this->_obs_double) {
                    this->_obs_double = std::unique_ptr<std::vector<Observable<double>>>(new std::vector<Observable<double>>());
                }
                Observable<double> superfluid_stiffness;
                superfluid_stiffness.set_observable_name("superfluid_stiffness");
                superfluid_stiffness.add_method(Methods::measure_superfluid_stiffness);
                this->_obs_double->push_back(superfluid_stiffness);
            }
        }


        // adding sign measurements manually to keep track of sign problem
        if (this->_obs_list_eqtime) {
            this->_sign_eqtime = std::unique_ptr<Observable<double>>(new Observable<double>);
            this->_sign_eqtime->set_observable_name("sign_eqtime");
            this->_sign_eqtime->add_method(Methods::measure_config_sign_eqtime);
        }
        if (this->_obs_list_dynamic) {
            this->_sign_dynamic = std::unique_ptr<Observable<double>>(new Observable<double>);
            this->_sign_dynamic->set_observable_name("sign_dynamic");
            this->_sign_dynamic->add_method(Methods::measure_config_sign_dynamic);
        }

        // adding new methods here
    }

} // namespace Measure
