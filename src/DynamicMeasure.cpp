#include "DynamicMeasure.h"
#include "Hubbard.h"

Measure::DynamicMeasure::DynamicMeasure(const int &nbin) {
    this->nbin = nbin;
}

void Measure::DynamicMeasure::resize(const int &_nbin) {
    this->nbin = _nbin;
}

void Measure::DynamicMeasure::initial(const Model::Hubbard &hubbard) {
    this->n_time_displaced = 0;
    this->sign.set_size_of_bin(this->nbin);
    this->superfluid_density.set_size_of_bin(this->nbin);

    this->matsubara_greens.reserve(hubbard.lt);
    this->density_of_states.reserve(hubbard.lt);
    for (int l = 0; l < hubbard.lt; ++l) {
        this->matsubara_greens.emplace_back(this->nbin);
        this->density_of_states.emplace_back(this->nbin);
    }
}

void Measure::DynamicMeasure::clear_temporary(const Model::Hubbard &hubbard) {
    this->n_time_displaced = 0;
    this->sign.clear_temporary();
    this->superfluid_density.clear_temporary();
    for (auto &green : this->matsubara_greens) {
        green.clear_temporary();
    }
    for (auto &dos : this->density_of_states) {
        dos.clear_temporary();
    }
}

void Measure::DynamicMeasure::time_displaced_measure(const Model::Hubbard &hubbard) {
    this->sign.tmp_value() += hubbard.config_sign;
    for (int l = 0; l < hubbard.lt; ++l) {
        this->measure_matsubara_greens(l, hubbard);
        this->measure_density_of_states(l, hubbard);
    }
    this->measure_superfluid_density(hubbard);
    this->n_time_displaced++;
}

void Measure::DynamicMeasure::normalize_stats(const Model::Hubbard &hubbard) {
    this->sign.tmp_value() /= this->n_time_displaced;
    this->superfluid_density.tmp_value() /= this->n_time_displaced * this->sign.tmp_value();
    for (auto &green : this->matsubara_greens) {
        green.tmp_value() /= this->n_time_displaced * this->sign.tmp_value();
    }
    for (auto &dos : this->density_of_states) {
        dos.tmp_value() /= this->n_time_displaced * this->sign.tmp_value();
    }
}

void Measure::DynamicMeasure::write_stats_to_bins(const int &bin, const Model::Hubbard &hubbard) {
    this->sign.bin_data()[bin] = this->sign.tmp_value();
    this->superfluid_density.bin_data()[bin] = this->superfluid_density.tmp_value();
    for (auto &green : this->matsubara_greens) {
        green.bin_data()[bin] = green.tmp_value();
    }
    for (auto &dos : this->density_of_states) {
        dos.bin_data()[bin] = dos.tmp_value();
    }
}

void Measure::DynamicMeasure::measure_matsubara_greens(const int &t, const Model::Hubbard &hubbard) {
    // factor 1/2 comes from two degenerate spin states
    const Eigen::MatrixXd gt0 = 0.5 * (hubbard.vec_green_t0_up[t] + hubbard.vec_green_t0_up[t]);
    double tmpFourier = 0.0;

    for (int xi = 0; xi < hubbard.ll; ++xi) {
        for (int yi = 0; yi < hubbard.ll; ++yi) {
            for (int xj = 0; xj < hubbard.ll; ++xj) {
                for (int yj = 0; yj < hubbard.ll; ++yj) {
                    const int i = xi + hubbard.ll * yi;
                    const int j = xj + hubbard.ll * yj;
                    const Eigen::VectorXd r = (Eigen::VectorXd(2) << (xi-xj), (yi-yj)).finished();
                    this->matsubara_greens[t].tmp_value() += hubbard.config_sign * cos(-r.dot(this->q)) * gt0(j, i) / hubbard.ls;
                }
            }
        }
    }
}

void Measure::DynamicMeasure::measure_density_of_states(const int &t, const Model::Hubbard &hubbard) {
    // spin degenerate model
    const Eigen::MatrixXd gt0 = 0.5 * (hubbard.vec_green_t0_up[t] + hubbard.vec_green_t0_up[t]);
    this->density_of_states[t].tmp_value() += hubbard.config_sign * gt0.trace() / hubbard.ls;
}

void Measure::DynamicMeasure::measure_superfluid_density(const Model::Hubbard &hubbard) {
    // momentum qx and qy
    const Eigen::VectorXd qx = (Eigen::VectorXd(2) << 2 * M_PI / hubbard.ll, 0).finished();
    const Eigen::VectorXd qy = (Eigen::VectorXd(2) << 0, 2 * M_PI / hubbard.ll).finished();

    // fourier transformation in time-energy space
    double tmp_rho_s = 0.0;
    for (int l = 0; l < hubbard.lt; ++l) {
        const Eigen::MatrixXd gt0_up = hubbard.vec_green_t0_up[l];
        const Eigen::MatrixXd g0t_up = hubbard.vec_green_0t_up[l];
        const Eigen::MatrixXd gtt_up = hubbard.vec_green_tt_up[l];
        const Eigen::MatrixXd g00_up = hubbard.vec_green_tt_up[hubbard.lt - 1];
        const Eigen::MatrixXd gt0_dn = hubbard.vec_green_t0_dn[l];
        const Eigen::MatrixXd g0t_dn = hubbard.vec_green_0t_dn[l];
        const Eigen::MatrixXd gtt_dn = hubbard.vec_green_tt_dn[l];
        const Eigen::MatrixXd g00_dn = hubbard.vec_green_tt_dn[hubbard.lt - 1];

        for (int xi = 0; xi < hubbard.ll; ++xi) {
            for (int yi = 0; yi < hubbard.ll; ++yi) {
                const int i = xi + hubbard.ll * yi;
                const int ipx = (xi + 1) % hubbard.ll + hubbard.ll * yi;

                for (int xj = 0; xj < hubbard.ll; ++xj) {
                    for (int yj = 0; yj < hubbard.ll; ++yj) {
                        /* for a given site l and time-slice tau
                         * the current-current correlation Jx-Jx: \Gamma_xx (l, \tau) = < jx(l, \tau) * jx(0, 0) > */
                        const int j = xj + hubbard.ll * yj;
                        const int jpx = (xj + 1) % hubbard.ll + hubbard.ll * yj;
                        const Eigen::VectorXd r = (Eigen::VectorXd(2) << (xi - xj), (yi - yj)).finished();
                        const double factor = hubbard.config_sign * cos(r.dot(qx)) - cos(r.dot(qy));

                        const double delta_tau_i_j = (l == 0 && i == j)? 1 : 0;
                        const double delta_tau_i_jpx = (l == 0 && i == jpx)? 1 : 0;
                        const double delta_tau_ipx_j = (l == 0 && ipx == j)? 1 : 0;

                        tmp_rho_s += hubbard.t * hubbard.t / hubbard.ls * factor * (
                                - ( gtt_up(ipx, i) - gtt_up(i, ipx) + gtt_dn(ipx, i) - gtt_dn(i, ipx) ) *
                                  ( g00_up(jpx, j) - g00_up(j, jpx) + g00_dn(jpx, j) - g00_dn(j, jpx) )

                                + gt0_up(ipx, jpx) * ( delta_tau_i_j - g0t_up(j, i) ) + gt0_dn(ipx, jpx) * ( delta_tau_i_j - g0t_dn(j, i))
                                - gt0_up(ipx, j) * ( delta_tau_i_jpx - g0t_up(jpx, i) ) + gt0_dn(ipx, j) * ( delta_tau_i_jpx - g0t_dn(jpx, i))
                                - gt0_up(i, jpx) * ( delta_tau_ipx_j - g0t_up(j, ipx) ) + gt0_dn(i, jpx) * ( delta_tau_ipx_j - g0t_dn(j, ipx))
                                + gt0_up(i, j) * ( delta_tau_i_j - g0t_up(jpx, ipx) ) + gt0_dn(i, j) * ( delta_tau_i_j - g0t_dn(jpx, ipx))
                                );
                    }
                }
            }
        }
    }
    this->superfluid_density.tmp_value() += 0.25 * tmp_rho_s;
}

void Measure::DynamicMeasure::analyse_stats(const Model::Hubbard &hubbard) {
    this->sign.analyse();
    this->superfluid_density.analyse();
    for (auto &green : this->matsubara_greens) {
        green.analyse();
    }
    for (auto &dos : this->density_of_states) {
        dos.analyse();
    }
}
