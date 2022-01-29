#include "measure.h"
#include "hubbard.h"

// Todo: passing vectors of pointers to avoid unnecessary copies

namespace Measure {

    int Measure::nbin() const { return this->_nbin; }
    bool Measure::is_eqtime_measure() const { return this->_container.is_eqtime_measure(); }
    bool Measure::is_dynamic_measure() const { return this->_container.is_dynamic_measure(); }

    std::vector<Observable<double>> Measure::obs_list_double() const { return this->_container.obs_list_double(); }
    std::vector<Observable<Eigen::VectorXd>> Measure::obs_list_vector() const { return this->_container.obs_list_vector(); }
    std::vector<Observable<Eigen::MatrixXd>> Measure::obs_list_matrix() const { return this->_container.obs_list_matrix(); }
    
    Observable<double> Measure::sign_eqtime() const { 
        return ((this->is_eqtime_measure())? *this->_container._sign_eqtime : Observable<double>()); 
    }
    Observable<double> Measure::sign_dynamic() const { 
        return ((this->is_dynamic_measure())? *this->_container._sign_dynamic : Observable<double>()); 
    }

    bool Measure::find(std::string obs_name) const {
        if (this->_container.is_eqtime_obs(obs_name) || this->_container.is_dynamic_obs(obs_name)) {
            return true;
        }
        else { return false; }
        // return ((this->_container.is_eqtime_obs(obs_name) || this->_container.is_dynamic_obs(obs_name))? true : false);
    }

    Observable<double> Measure::find_double_obs(std::string obs_name) const {
        for (auto obs : this->obs_list_double()) {
            if (obs_name == obs.name()) {
                return obs;
            }
        }
        return Observable<double>();
    }

    Observable<Eigen::VectorXd> Measure::find_vector_obs(std::string obs_name) const {
        for (auto obs : this->obs_list_vector()) {
            if (obs_name == obs.name()) {
                return obs;
            }
        }
        return Observable<Eigen::VectorXd>();
    }

    Observable<Eigen::MatrixXd> Measure::find_matrix_obs(std::string obs_name) const {
        for (auto obs : this->obs_list_matrix()) {
            if (obs_name == obs.name()) {
                return obs;
            }
        }
        return Observable<Eigen::MatrixXd>();
    }

    void Measure::set_size_of_bin(const int &nbin) {
        this->_nbin = nbin;
    }

    void Measure::set_observable_list(const std::vector<std::string> &obs_list) {
        this->_obs_list = obs_list;
    }

    void Measure::set_lattice_momentum(const Eigen::VectorXd &q) {
        assert( q.size() == 2 );
        this->q = q;
    }

    void Measure::initial(const Model::Hubbard &hubbard) {
        // initialize container
        this->_container.read_list(this->_obs_list);
        this->_container.filter();
        this->_container.allocate();

        // allocate memory
        if (this->_container._obs_double) {
            for (auto &obs_double : *(this->_container._obs_double)) {
                obs_double.set_zero_element(0.0);
                obs_double.set_size_of_bin(this->_nbin);
                obs_double.allocate();
            }
        }
        
        if (this->_container._obs_vector) {
            for (auto &obs_vector : *(this->_container._obs_vector)) {
                obs_vector.set_zero_element(Eigen::VectorXd::Zero(hubbard.lt));
                obs_vector.set_size_of_bin(this->_nbin);
                obs_vector.allocate();
            }
        }
        
        if (this->_container._obs_matrix) {
            for (auto &obs_matrix : *(this->_container._obs_matrix)) {
                obs_matrix.set_zero_element(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
                obs_matrix.set_size_of_bin(this->_nbin);
                obs_matrix.allocate();
            }
        }

        // sign measurements
        if (this->_container._sign_eqtime) {
            this->_container._sign_eqtime->set_zero_element(0.0);
            this->_container._sign_eqtime->set_size_of_bin(this->_nbin);
            this->_container._sign_eqtime->allocate();
        }
        if (this->_container._sign_dynamic) {
            this->_container._sign_dynamic->set_zero_element(0.0);
            this->_container._sign_dynamic->set_size_of_bin(this->_nbin);
            this->_container._sign_dynamic->allocate();
        }
        
    }

    void Measure::clear_temporary() {
        if (this->_container._obs_double) {
            for (auto &obs_double : *(this->_container._obs_double)) {
                obs_double.clear_temporary();
            }
        }
        if (this->_container._obs_vector) {
            for (auto &obs_vector : *(this->_container._obs_vector)) {
                obs_vector.clear_temporary();
            }
        }
        if (this->_container._obs_matrix) {
            for (auto &obs_matrix : *(this->_container._obs_matrix)) {
                obs_matrix.clear_temporary();
            }
        }
        if (this->_container._sign_eqtime) { this->_container._sign_eqtime->clear_temporary(); }
        if (this->_container._sign_dynamic) { this->_container._sign_dynamic->clear_temporary(); }
    }

    void Measure::eqtime_measure(const Model::Hubbard &hubbard) {
        if (this->_container._sign_eqtime) {
            this->_container._sign_eqtime->measure(*this, hubbard);
        }

        // Todo: improve efficiency here
        if (this->_container._obs_double) {
            for (auto &obs_double : *(this->_container._obs_double)) {
                if (this->_container.is_eqtime_obs(obs_double.name())) {
                    obs_double.measure(*this, hubbard);
                }
            }
        }
        if (this->_container._obs_vector) {
            for (auto &obs_vector : *(this->_container._obs_vector)) {
                if (this->_container.is_eqtime_obs(obs_vector.name())) {
                    obs_vector.measure(*this, hubbard);
                }
            }
        }
        if (this->_container._obs_matrix) {
            for (auto &obs_matrix : *(this->_container._obs_matrix)) {
                if (this->_container.is_eqtime_obs(obs_matrix.name())) {
                    obs_matrix.measure(*this, hubbard);
                }
            }
        }
    }

    void Measure::dynamic_measure(const Model::Hubbard &hubbard) {
        if (this->_container._sign_dynamic) {
            this->_container._sign_dynamic->measure(*this, hubbard);
        }

        // Todo: improve efficiency here
        if (this->_container._obs_double) {
            for (auto &obs_double : *(this->_container._obs_double)) {
                if (this->_container.is_dynamic_obs(obs_double.name())) {
                    obs_double.measure(*this, hubbard);
                }
            }
        }
        if (this->_container._obs_vector) {
            for (auto &obs_vector : *(this->_container._obs_vector)) {
                if (this->_container.is_dynamic_obs(obs_vector.name())) {
                    obs_vector.measure(*this, hubbard);
                }
            }
        }
        if (this->_container._obs_matrix) {
            for (auto &obs_matrix : *(this->_container._obs_matrix)) {
                if (this->_container.is_dynamic_obs(obs_matrix.name())) {
                    obs_matrix.measure(*this, hubbard);
                }
            }
        }
    }

    void Measure::normalize_stats() {
        if (this->_container._sign_eqtime) {
            this->_container._sign_eqtime->tmp_value() /= this->_container._sign_eqtime->counts();
        }
        if (this->_container._sign_dynamic) {
            this->_container._sign_dynamic->tmp_value() /= this->_container._sign_dynamic->counts();
        }

        if (this->_container._obs_double) {
            for (auto &obs_double : *(this->_container._obs_double)) {
                if ( this->_container.is_eqtime_obs(obs_double.name()) ) {
                    obs_double.tmp_value() /= obs_double.counts() * this->_container._sign_eqtime->tmp_value();
                }
                if ( this->_container.is_dynamic_obs(obs_double.name()) ) {
                    obs_double.tmp_value() /= obs_double.counts() * this->_container._sign_dynamic->tmp_value();
                }
            }
        }

        if (this->_container._obs_vector) {
            for (auto &obs_vector : *(this->_container._obs_vector)) {
                if ( this->_container.is_eqtime_obs(obs_vector.name()) ) {
                    obs_vector.tmp_value() /= obs_vector.counts() * this->_container._sign_eqtime->tmp_value();
                }
                if ( this->_container.is_dynamic_obs(obs_vector.name()) ) {
                    obs_vector.tmp_value() /= obs_vector.counts() * this->_container._sign_dynamic->tmp_value();
                }
            }
        }

        if (this->_container._obs_matrix) {
            for (auto &obs_matrix : *(this->_container._obs_matrix)) {
                if ( this->_container.is_eqtime_obs(obs_matrix.name()) ) {
                    obs_matrix.tmp_value() /= obs_matrix.counts() * this->_container._sign_eqtime->tmp_value();
                }
                if ( this->_container.is_dynamic_obs(obs_matrix.name()) ) {
                    obs_matrix.tmp_value() /= obs_matrix.counts() * this->_container._sign_dynamic->tmp_value();
                }
            }
        }
    }

    void Measure::write_stats_to_bins(int bin) {
        if (this->_container._sign_eqtime) {
            this->_container._sign_eqtime->bin_data()[bin] = this->_container._sign_eqtime->tmp_value();
        }
        if (this->_container._sign_dynamic) {
            this->_container._sign_dynamic->bin_data()[bin] = this->_container._sign_dynamic->tmp_value();
        }

        if (this->_container._obs_double) {
            for (auto &obs_double : *(this->_container._obs_double)) {
                obs_double.bin_data()[bin] = obs_double.tmp_value();
            }
        }

        if (this->_container._obs_vector) {
            for (auto &obs_vector : *(this->_container._obs_vector)) {
                obs_vector.bin_data()[bin] = obs_vector.tmp_value();
            }
        }

        if (this->_container._obs_matrix) {
            for (auto &obs_matrix : *(this->_container._obs_matrix)) {
                obs_matrix.bin_data()[bin] = obs_matrix.tmp_value();
            }
        }
    }

    void Measure::analyse_stats() {
        if (this->_container._sign_eqtime) {
            this->_container._sign_eqtime->analyse();
        }
        if (this->_container._sign_dynamic) {
            this->_container._sign_dynamic->analyse();
        }

        if (this->_container._obs_double) {
            for (auto &obs_double : *(this->_container._obs_double)) {
                obs_double.analyse();
            }
        }

        if (this->_container._obs_vector) {
            for (auto &obs_vector : *(this->_container._obs_vector)) {
                obs_vector.analyse();
            }
        }

        if (this->_container._obs_matrix) {
            for (auto &obs_matrix : *(this->_container._obs_matrix)) {
                obs_matrix.analyse();
            }
        }
    }

} // namespace Measure