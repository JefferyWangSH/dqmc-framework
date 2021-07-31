#include "eqtimeMeasure.h"
#include "hubbard.h"

void measure::eqtimeMeasure::resize(const int &nbin) {
    this->nbin = nbin;
}

void measure::eqtimeMeasure::initial() {
    obs_bin_eqtime["double_occupancy"].reserve(nbin);
    obs_bin_eqtime["kinetic_energy"].reserve(nbin);
    obs_bin_eqtime["structure_factor"].reserve(nbin);
    obs_bin_eqtime["momentum_distribution"].reserve(nbin);
    obs_bin_eqtime["local_spin_correlation"].reserve(nbin);
    obs_bin_eqtime["average_sign"].reserve(nbin);
}

void measure::eqtimeMeasure::clear() {
    n_equal_time = 0;
    double_occupancy = 0.0;
    kinetic_energy = 0.0;
    structure_factor = 0.0;
    momentum_distribution = 0.0;
    local_spin_correlation = 0.0;
    average_sign = 0.0;
}

void measure::eqtimeMeasure::meas_Double_Occu(const Hubbard &hubbard, const int &t) {
    assert( t >= 0 && t < hubbard.lt);
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];

    for (int i = 0; i < hubbard.ls; ++i) {
        const double tmp_double_occupancy = (1 - gu(i,i)) * (1 - gd(i,i));
        double_occupancy += hubbard.config_sign * tmp_double_occupancy;
    }
}

void measure::eqtimeMeasure::meas_Kinetic_Energy(const Hubbard &hubbard, const int &t) {
    assert( t >= 0 && t < hubbard.lt);
    const int ll = hubbard.ll;
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];

    for (int x = 0; x < ll; ++x) {
        for (int y = 0; y < ll; ++y) {
            const double tmp_kinetic_energy = 2 * hubbard.t * (gu(x + ll*y, ((x+1)%ll) + ll*y) + gu(x + ll*y, x + ll*((y+1)%ll)))
                                   + 2 * hubbard.t * (gd(x + ll*y, ((x+1)%ll) + ll*y) + gd(x + ll*y, x + ll*((y+1)%ll)));
            kinetic_energy += hubbard.config_sign * tmp_kinetic_energy;
        }
    }
}

void measure::eqtimeMeasure::meas_Momentum_Dist(const Hubbard &hubbard, const int &t, const Eigen::VectorXd &p) {
    assert( t >= 0 && t < hubbard.lt);
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
    momentum_distribution += hubbard.config_sign * (1 - 0.5 * tmp_fourier / hubbard.ls);
}

void measure::eqtimeMeasure::meas_local_Spin_Corr(const Hubbard &hubbard, const int &t) {
    assert( t >= 0 && t < hubbard.lt);
    const int ls = hubbard.ls;
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];
    double  tmp_onsite_correlation = 0.0;

    for (int i = 0; i < ls; ++i) {
        tmp_onsite_correlation += gu(i, i) + gd(i, i) - 2 * gu(i, i) * gd(i, i);
    }
    local_spin_correlation += hubbard.config_sign * tmp_onsite_correlation / ls;
}

void measure::eqtimeMeasure::meas_Struct_Factor(const Hubbard &hubbard, const int &t, const Eigen::VectorXd &p) {
    assert( t >= 0 && t < hubbard.lt);
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
                    structure_factor += hubbard.config_sign * tmp_structure_factor;
                }
            }
        }
    }
}

void measure::eqtimeMeasure::measure_equal_time(const Hubbard &hubbard) {
    for (int t = 0; t < hubbard.lt; ++t) {
        meas_Double_Occu(hubbard, t);
        meas_Kinetic_Energy(hubbard, t);
        meas_Struct_Factor(hubbard, t, q);
        meas_Momentum_Dist(hubbard, t, q);
        meas_local_Spin_Corr(hubbard, t);
    }
    average_sign += hubbard.config_sign;
    n_equal_time++;
}

void measure::eqtimeMeasure::normalizeStats(const Hubbard &hubbard) {
    average_sign /= n_equal_time;
    double_occupancy /= hubbard.ls * hubbard.lt * n_equal_time * average_sign;
    kinetic_energy /= hubbard.ls * hubbard.lt * n_equal_time * average_sign;
    structure_factor /= hubbard.ls * hubbard.ls * hubbard.lt * n_equal_time * average_sign;
    momentum_distribution /= hubbard.lt * n_equal_time * average_sign;
    local_spin_correlation /= hubbard.lt * n_equal_time * average_sign;
}

void measure::eqtimeMeasure::write_Stats_to_bins(int bin) {
    obs_bin_eqtime["double_occupancy"][bin] = double_occupancy;
    obs_bin_eqtime["kinetic_energy"][bin] = kinetic_energy;
    obs_bin_eqtime["structure_factor"][bin] = structure_factor;
    obs_bin_eqtime["momentum_distribution"][bin] = momentum_distribution;
    obs_bin_eqtime["local_spin_correlation"][bin] = local_spin_correlation;
    obs_bin_eqtime["average_sign"][bin] = average_sign;
}

void measure::eqtimeMeasure::analyse_equal_time_Stats(const std::string &obs) {
    assert(obs_bin_eqtime.count(obs) == 1);

    // clear data of previous statistics
    obs_mean_eqtime[obs] = 0;
    obs_err_eqtime[obs] = 0;

    for (int bin = 0; bin < nbin; ++bin) {
        obs_mean_eqtime[obs] += obs_bin_eqtime[obs][bin];
        obs_err_eqtime[obs] += pow(obs_bin_eqtime[obs][bin], 2);
    }

    obs_mean_eqtime[obs] /= nbin;
    obs_err_eqtime[obs] /= nbin;
    obs_err_eqtime[obs] = pow(obs_err_eqtime[obs] - pow(obs_mean_eqtime[obs], 2), 0.5) / pow(nbin - 1, 0.5);
}

void measure::eqtimeMeasure::analyseStats() {
    analyse_equal_time_Stats("double_occupancy");
    analyse_equal_time_Stats("kinetic_energy");
    analyse_equal_time_Stats("structure_factor");
    analyse_equal_time_Stats("momentum_distribution");
    analyse_equal_time_Stats("local_spin_correlation");
    analyse_equal_time_Stats("average_sign");
}
