#include "EqtimeMeasure.h"
#include "Hubbard.h"

Measure::EqtimeMeasure::EqtimeMeasure(const int &nbin) {
    this->nbin = nbin;
}

void Measure::EqtimeMeasure::resize(const int &_nbin) {
    this->nbin = _nbin;
}

void Measure::EqtimeMeasure::initial() {
    this->obs_bin_eqtime["double_occupancy"].reserve(nbin);
    this->obs_bin_eqtime["kinetic_energy"].reserve(nbin);
    this->obs_bin_eqtime["structure_factor"].reserve(nbin);
    this->obs_bin_eqtime["momentum_distribution"].reserve(nbin);
    this->obs_bin_eqtime["local_spin_correlation"].reserve(nbin);
    this->obs_bin_eqtime["average_sign"].reserve(nbin);
}

void Measure::EqtimeMeasure::clear() {
    this->n_equal_time = 0;
    this->double_occupancy = 0.0;
    this->kinetic_energy = 0.0;
    this->structure_factor = 0.0;
    this->momentum_distribution = 0.0;
    this->local_spin_correlation = 0.0;
    this->average_sign = 0.0;
}

void Measure::EqtimeMeasure::meas_Double_Occu(const Model::Hubbard &hubbard, const int &t) {
    assert( t >= 0 && t < hubbard.lt );
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];

    for (int i = 0; i < hubbard.ls; ++i) {
        const double tmp_double_occupancy = (1 - gu(i,i)) * (1 - gd(i,i));
        this->double_occupancy += hubbard.config_sign * tmp_double_occupancy;
    }
}

void Measure::EqtimeMeasure::meas_Kinetic_Energy(const Model::Hubbard &hubbard, const int &t) {
    assert( t >= 0 && t < hubbard.lt );
    const int ll = hubbard.ll;
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];

    for (int x = 0; x < ll; ++x) {
        for (int y = 0; y < ll; ++y) {
            const double tmp_kinetic_energy = 2 * hubbard.t * (gu(x + ll*y, ((x+1)%ll) + ll*y) + gu(x + ll*y, x + ll*((y+1)%ll)))
                                   + 2 * hubbard.t * (gd(x + ll*y, ((x+1)%ll) + ll*y) + gd(x + ll*y, x + ll*((y+1)%ll)));
            this->kinetic_energy += hubbard.config_sign * tmp_kinetic_energy;
        }
    }
}

void Measure::EqtimeMeasure::meas_Momentum_Dist(const Model::Hubbard &hubbard, const int &t, const Eigen::VectorXd &p) {
    assert( t >= 0 && t < hubbard.lt );
    const int ll = hubbard.ll;
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];
    double tmp_fourier = 0.0;

    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            for (int xj = 0; xj < ll; ++xj) {
                for (int yj = 0; yj < ll; ++yj) {
                    const int i = xi + ll * yi;
                    const int j = xj + ll * yj;
                    const Eigen::VectorXd r = (Eigen::VectorXd(2) << (xi-xj), (yi-yj)).finished();
                    tmp_fourier += cos(-r.dot(p)) * (gu(j, i) + gd(j, i));
                }
            }
        }
    }
    this->momentum_distribution += hubbard.config_sign * (1 - 0.5 * tmp_fourier / hubbard.ls);
}

void Measure::EqtimeMeasure::meas_local_Spin_Corr(const Model::Hubbard &hubbard, const int &t) {
    assert( t >= 0 && t < hubbard.lt );
    const int ls = hubbard.ls;
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];
    double tmp_onsite_correlation = 0.0;

    for (int i = 0; i < ls; ++i) {
        tmp_onsite_correlation += gu(i, i) + gd(i, i) - 2 * gu(i, i) * gd(i, i);
    }
    this->local_spin_correlation += hubbard.config_sign * tmp_onsite_correlation / ls;
}

void Measure::EqtimeMeasure::meas_Struct_Factor(const Model::Hubbard &hubbard, const int &t, const Eigen::VectorXd &p) {
    assert( t >= 0 && t < hubbard.lt );
    const int ll = hubbard.ll;
    const int ls = hubbard.ls;
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];

    /**  gu(i,j) = < c_i c^+_j >
     *  guc(i,j) = < c^+_i c_j > */
    Eigen::MatrixXd guc = Eigen::MatrixXd::Identity(ls, ls);
    Eigen::MatrixXd gdc = Eigen::MatrixXd::Identity(ls, ls);

    // get guc and gdc
    for (int i = 0; i < ls; ++i) {
        for (int j = 0; j < ls; ++j) {
            guc(j, i) = - gu(i, j);
            gdc(j, i) = - gd(i, j);
        }
        guc(i, i)++;
        gdc(i, i)++;
    }

    // loop for site i, j
    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            for (int xj = 0; xj < ll; ++xj) {
                for (int yj = 0; yj < ll; ++yj) {
                    const int i = xi + ll * yi;
                    const int j = xj + ll * yj;
                    const Eigen::VectorXd r = (Eigen::VectorXd(2) << (xi-xj), (yi-yj)).finished();
                    const double factor = cos(-r.dot(p));
                    /** factor 0.25 comes from spin 1/2 */
                    const double tmp_structure_factor = 0.25 * factor * (
                            + guc(i, i) * guc(j, j) + guc(i, j) * gu(i, j)
                            + gdc(i, i) * gdc(j, j) + gdc(i, j) * gd(i, j)
                            - gdc(i, i) * guc(j, j) - guc(i, i) * gdc(j, j)
                            );
                    this->structure_factor += hubbard.config_sign * tmp_structure_factor;
                }
            }
        }
    }
}

void Measure::EqtimeMeasure::measure_equal_time(const Model::Hubbard &hubbard) {
    for (int t = 0; t < hubbard.lt; ++t) {
        this->meas_Double_Occu(hubbard, t);
        this->meas_Kinetic_Energy(hubbard, t);
        this->meas_Struct_Factor(hubbard, t, this->q);
        this->meas_Momentum_Dist(hubbard, t, this->q);
        this->meas_local_Spin_Corr(hubbard, t);
    }
    this->average_sign += hubbard.config_sign;
    this->n_equal_time++;
}

void Measure::EqtimeMeasure::normalizeStats(const Model::Hubbard &hubbard) {
    this->average_sign /= this->n_equal_time;
    this->double_occupancy /= hubbard.ls * hubbard.lt * this->n_equal_time * this->average_sign;
    this->kinetic_energy /= hubbard.ls * hubbard.lt * this->n_equal_time * this->average_sign;
    this->structure_factor /= hubbard.ls * hubbard.ls * hubbard.lt * this->n_equal_time * this->average_sign;
    this->momentum_distribution /= hubbard.lt * this->n_equal_time * this->average_sign;
    this->local_spin_correlation /= hubbard.lt * this->n_equal_time * this->average_sign;
}

void Measure::EqtimeMeasure::write_Stats_to_bins(int bin) {
    this->obs_bin_eqtime["double_occupancy"][bin] = this->double_occupancy;
    this->obs_bin_eqtime["kinetic_energy"][bin] = this->kinetic_energy;
    this->obs_bin_eqtime["structure_factor"][bin] = this->structure_factor;
    this->obs_bin_eqtime["momentum_distribution"][bin] = this->momentum_distribution;
    this->obs_bin_eqtime["local_spin_correlation"][bin] = this->local_spin_correlation;
    this->obs_bin_eqtime["average_sign"][bin] = this->average_sign;
}

void Measure::EqtimeMeasure::analyse_equal_time_Stats(const std::string &obs) {
    assert( obs_bin_eqtime.count(obs) == 1 );

    // clear data of previous statistics
    this->obs_mean_eqtime[obs] = 0;
    this->obs_err_eqtime[obs] = 0;

    for (int bin = 0; bin < this->nbin; ++bin) {
        this->obs_mean_eqtime[obs] += this->obs_bin_eqtime[obs][bin];
        this->obs_err_eqtime[obs] += pow(this->obs_bin_eqtime[obs][bin], 2);
    }

    this->obs_mean_eqtime[obs] /= this->nbin;
    this->obs_err_eqtime[obs] /= this->nbin;
    this->obs_err_eqtime[obs] = pow(this->obs_err_eqtime[obs] - pow(this->obs_mean_eqtime[obs], 2), 0.5) / pow(this->nbin - 1, 0.5);
}

void Measure::EqtimeMeasure::analyseStats() {
    this->analyse_equal_time_Stats("double_occupancy");
    this->analyse_equal_time_Stats("kinetic_energy");
    this->analyse_equal_time_Stats("structure_factor");
    this->analyse_equal_time_Stats("momentum_distribution");
    this->analyse_equal_time_Stats("local_spin_correlation");
    this->analyse_equal_time_Stats("average_sign");
}
