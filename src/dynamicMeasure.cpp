#include "dynamicMeasure.h"
#include "hubbard.h"

void measure::dynamicMeasure::resize(const int &nbin) {
    this->nbin = nbin;
}

void measure::dynamicMeasure::initial(const Hubbard &hubbard) {
    // clear data of previous simulation
    clear();

    obs_bin_g_kt.reserve(nbin);
    obs_bin_rho_s.reserve(nbin);
    obs_bin_sign.reserve(nbin);

    obs_bin_gt0_up.reserve(nbin);
    obs_bin_g0t_up.reserve(nbin);
    obs_bin_gtt_up.reserve(nbin);
    obs_bin_gt0_dn.reserve(nbin);
    obs_bin_g0t_dn.reserve(nbin);
    obs_bin_gtt_dn.reserve(nbin);

    for (int bin = 0; bin < nbin; ++bin) {
        obs_bin_g_kt.emplace_back(hubbard.lt, 0.0);
        obs_bin_rho_s.emplace_back(0.0);
        obs_bin_sign.emplace_back(0.0);

        obs_bin_gt0_up.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        obs_bin_g0t_up.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        obs_bin_gtt_up.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        obs_bin_gt0_dn.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        obs_bin_g0t_dn.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        obs_bin_gtt_dn.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
    }

    obs_mean_rho_s = 0.0;
    obs_err_rho_s = 0.0;
    
    obs_mean_sign = 0.0;
    obs_err_sign = 0.0;
    tmp_sign = 0.0;

    obs_mean_g_kt.reserve(hubbard.lt);
    obs_err_g_kt.reserve(hubbard.lt);

    tmp_gt0_tau_up.reserve(hubbard.lt);
    tmp_g0t_tau_up.reserve(hubbard.lt);
    tmp_gtt_tau_up.reserve(hubbard.lt);
    tmp_gt0_tau_dn.reserve(hubbard.lt);
    tmp_g0t_tau_dn.reserve(hubbard.lt);
    tmp_gtt_tau_dn.reserve(hubbard.lt);

    for (int l = 0; l < hubbard.lt; ++l) {
        obs_mean_g_kt.emplace_back(0.0);
        obs_err_g_kt.emplace_back(0.0);

        tmp_gt0_tau_up.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        tmp_g0t_tau_up.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        tmp_gtt_tau_up.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        tmp_gt0_tau_dn.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        tmp_g0t_tau_dn.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        tmp_gtt_tau_dn.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
    }

    // free redundant memory
    obs_bin_g_kt.shrink_to_fit();
    obs_mean_g_kt.shrink_to_fit();
    obs_err_g_kt.shrink_to_fit();
    obs_bin_rho_s.shrink_to_fit();
    obs_bin_sign.shrink_to_fit();

    obs_bin_gt0_up.shrink_to_fit();
    obs_bin_g0t_up.shrink_to_fit();
    obs_bin_gtt_up.shrink_to_fit();
    obs_bin_gt0_dn.shrink_to_fit();
    obs_bin_g0t_dn.shrink_to_fit();
    obs_bin_gtt_dn.shrink_to_fit();

    tmp_gt0_tau_up.shrink_to_fit();
    tmp_g0t_tau_up.shrink_to_fit();
    tmp_gtt_tau_up.shrink_to_fit();
    tmp_gt0_tau_dn.shrink_to_fit();
    tmp_g0t_tau_dn.shrink_to_fit();
    tmp_gtt_tau_dn.shrink_to_fit();

    assert(obs_bin_gt0_up.size() == nbin);
    assert(obs_bin_g0t_up.size() == nbin);
    assert(obs_bin_gtt_up.size() == nbin);
    assert(obs_bin_gt0_dn.size() == nbin);
    assert(obs_bin_g0t_dn.size() == nbin);
    assert(obs_bin_gtt_dn.size() == nbin);

    for (int bin = 0; bin < nbin; ++bin) {
        obs_bin_gt0_up.shrink_to_fit();
        obs_bin_g0t_up.shrink_to_fit();
        obs_bin_gtt_up.shrink_to_fit();
        obs_bin_gt0_dn.shrink_to_fit();
        obs_bin_g0t_dn.shrink_to_fit();
        obs_bin_gtt_dn.shrink_to_fit();
    }
}

void measure::dynamicMeasure::clear() {
    obs_bin_g_kt.clear();
    obs_bin_rho_s.clear();
    obs_bin_sign.clear();

    obs_bin_gt0_up.clear();
    obs_bin_g0t_up.clear();
    obs_bin_gtt_up.clear();
    obs_bin_gt0_dn.clear();
    obs_bin_g0t_dn.clear();
    obs_bin_gtt_dn.clear();

    obs_mean_g_kt.clear();
    obs_err_g_kt.clear();

    obs_mean_rho_s = 0.0;
    obs_err_rho_s = 0.0;
    obs_mean_sign = 0.0;
    obs_err_sign = 0.0;

    tmp_gt0_tau_up.clear();
    tmp_g0t_tau_up.clear();
    tmp_gtt_tau_up.clear();
    tmp_gt0_tau_dn.clear();
    tmp_g0t_tau_dn.clear();
    tmp_gtt_tau_dn.clear();
}

void measure::dynamicMeasure::clear(const Hubbard &hubbard) {
    assert(tmp_gt0_tau_up.size() == hubbard.lt);
    assert(tmp_g0t_tau_up.size() == hubbard.lt);
    assert(tmp_gtt_tau_up.size() == hubbard.lt);
    assert(tmp_gt0_tau_dn.size() == hubbard.lt);
    assert(tmp_g0t_tau_dn.size() == hubbard.lt);
    assert(tmp_gtt_tau_dn.size() == hubbard.lt);

    n_time_displaced = 0;
    tmp_sign = 0.0;
    for (int l = 0; l < hubbard.lt; ++l) {
        tmp_gt0_tau_up[l].setZero(hubbard.ls, hubbard.ls);
        tmp_g0t_tau_up[l].setZero(hubbard.ls, hubbard.ls);
        tmp_gtt_tau_up[l].setZero(hubbard.ls, hubbard.ls);
        tmp_gt0_tau_dn[l].setZero(hubbard.ls, hubbard.ls);
        tmp_g0t_tau_dn[l].setZero(hubbard.ls, hubbard.ls);
        tmp_gtt_tau_dn[l].setZero(hubbard.ls, hubbard.ls);
    }
}

void measure::dynamicMeasure::measure_time_displaced(const Hubbard &hubbard) {
    n_time_displaced++;
    for (int l = 0; l < hubbard.lt; ++l) {
        tmp_gt0_tau_up[l] += hubbard.config_sign * hubbard.vec_green_t0_up[l];
        tmp_g0t_tau_up[l] += hubbard.config_sign * hubbard.vec_green_0t_up[l];
        tmp_gtt_tau_up[l] += hubbard.config_sign * hubbard.vec_green_tt_up[l];
        tmp_gt0_tau_dn[l] += hubbard.config_sign * hubbard.vec_green_t0_dn[l];
        tmp_g0t_tau_dn[l] += hubbard.config_sign * hubbard.vec_green_0t_dn[l];
        tmp_gtt_tau_dn[l] += hubbard.config_sign * hubbard.vec_green_tt_dn[l];
    }
    tmp_sign += hubbard.config_sign;
}

void measure::dynamicMeasure::normalizeStats(const Hubbard &hubbard) {
    for (int l = 0; l < hubbard.lt; ++l) {
        tmp_gt0_tau_up[l] /= n_time_displaced;
        tmp_g0t_tau_up[l] /= n_time_displaced;
        tmp_gtt_tau_up[l] /= n_time_displaced;
        tmp_gt0_tau_dn[l] /= n_time_displaced;
        tmp_g0t_tau_dn[l] /= n_time_displaced;
        tmp_gtt_tau_dn[l] /= n_time_displaced;
    }
    tmp_sign /= n_time_displaced;
}

void measure::dynamicMeasure::write_Stats_to_bins(const int &bin, const Hubbard &hubbard) {
    // TODO: CHECK IT OUT
    for (int l = 0; l < hubbard.lt; ++l) {
        obs_bin_gt0_up[bin][l] = tmp_gt0_tau_up[l] / tmp_sign;
        obs_bin_g0t_up[bin][l] = tmp_g0t_tau_up[l] / tmp_sign;
        obs_bin_gtt_up[bin][l] = tmp_gtt_tau_up[l] / tmp_sign;
        obs_bin_gt0_dn[bin][l] = tmp_gt0_tau_dn[l] / tmp_sign;
        obs_bin_g0t_dn[bin][l] = tmp_g0t_tau_dn[l] / tmp_sign;
        obs_bin_gtt_dn[bin][l] = tmp_gtt_tau_dn[l] / tmp_sign;
    }
    obs_bin_sign[bin] = tmp_sign;
}

void measure::dynamicMeasure::analyse_Dynamical_Corr(const int &bin, const Hubbard &hubbard) {
    for (int l = 0; l < hubbard.lt; ++l) {
        /* factor 0.5 comes from two spin states */
        const Eigen::MatrixXd gt0 = 0.5 * (obs_bin_gt0_up[bin][l] + obs_bin_gt0_dn[bin][l]);
        double tmpFourier = 0.0;

        for (int xi = 0; xi < hubbard.ll; ++xi) {
            for (int yi = 0; yi < hubbard.ll; ++yi) {
                for (int xj = 0; xj < hubbard.ll; ++xj) {
                    for (int yj = 0; yj < hubbard.ll; ++yj) {
                        const int i = xi + hubbard.ll * yi;
                        const int j = xj + hubbard.ll * yj;
                        const Eigen::VectorXd r = (Eigen::VectorXd(2) << (xi-xj), (yi-yj)).finished();
                        tmpFourier += cos(-r.dot(q)) * gt0(j, i) / hubbard.ls;
                    }
                }
            }
        }
        obs_bin_g_kt[bin][l] = tmpFourier;
    }
}

void measure::dynamicMeasure::analyse_Rho_S(const int &bin, const Hubbard &hubbard) {

    // momentum qx and qy
    const Eigen::VectorXd qx = (Eigen::VectorXd(2) << 2 * M_PI / hubbard.ll, 0).finished();
    const Eigen::VectorXd qy = (Eigen::VectorXd(2) << 0, 2 * M_PI / hubbard.ll).finished();

    double tmp_rho_s = 0.0;

    // fourier transformation in time-energy space
    for (int l = 0; l < hubbard.lt; ++l) {
        double tmp_fourier = 0.0;

        const Eigen::MatrixXd gt0_up = obs_bin_gt0_up[bin][l];
        const Eigen::MatrixXd g0t_up = obs_bin_g0t_up[bin][l];
        const Eigen::MatrixXd gtt_up = obs_bin_gtt_up[bin][l];
        const Eigen::MatrixXd g00_up = obs_bin_gtt_up[bin][hubbard.lt - 1];
        const Eigen::MatrixXd gt0_dn = obs_bin_gt0_dn[bin][l];
        const Eigen::MatrixXd g0t_dn = obs_bin_g0t_dn[bin][l];
        const Eigen::MatrixXd gtt_dn = obs_bin_gtt_dn[bin][l];
        const Eigen::MatrixXd g00_dn = obs_bin_gtt_dn[bin][hubbard.lt - 1];

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
                        const double factor = cos(r.dot(qx)) - cos(r.dot(qy));

                        const double delta_tau_i_j = (l == 0 && i == j)? 1 : 0;
                        const double delta_tau_i_jpx = (l == 0 && i == jpx)? 1 : 0;
                        const double delta_tau_ipx_j = (l == 0 && ipx == j)? 1 : 0;

                        tmp_fourier += hubbard.t * hubbard.t * factor * (
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
        tmp_fourier /= hubbard.ls * hubbard.ls;
        tmp_rho_s += tmp_fourier * hubbard.lt;      //FIXME: why this work?
    }

    obs_bin_rho_s[bin] = 0.25 * tmp_rho_s;
}

void measure::dynamicMeasure::analyse_timeDisplaced_Stats(const Hubbard &hubbard) {
    // analyse time-displaced statistics and record data.

    // clear previous statistics
    for (int l = 0; l < hubbard.lt; ++l) {
        obs_mean_g_kt[l] = 0.0;
        obs_err_g_kt[l] = 0.0;
    }
    obs_mean_rho_s = 0.0;
    obs_err_rho_s = 0.0;
    obs_mean_sign = 0.0;
    obs_err_sign = 0.0;

    // analyse statistics
    for (int bin = 0; bin < nbin; ++bin) {
        analyse_Dynamical_Corr(bin, hubbard);
        analyse_Rho_S(bin, hubbard);
    }

    // sum over bins and calculate mean and error
    for (int bin = 0; bin < nbin; ++bin) {
        for (int l = 0; l < hubbard.lt; ++l) {
            obs_mean_g_kt[l] += obs_bin_g_kt[bin][l];
            obs_err_g_kt[l] += obs_bin_g_kt[bin][l] * obs_bin_g_kt[bin][l];
        }
        obs_mean_rho_s += obs_bin_rho_s[bin];
        obs_err_rho_s += obs_bin_rho_s[bin] * obs_bin_rho_s[bin];
        obs_mean_sign += obs_bin_sign[bin];
        obs_err_sign += obs_bin_sign[bin] * obs_bin_sign[bin];
    }

    obs_mean_sign /= nbin;
    obs_err_sign /= nbin;
    obs_err_sign = pow(obs_err_sign - pow(obs_mean_sign, 2), 0.5) / pow(nbin - 1, 0.5);
    for (int l = 0; l < hubbard.lt; ++l) {
        obs_mean_g_kt[l] /= nbin;
        obs_err_g_kt[l] /= nbin;
        obs_err_g_kt[l] = pow(obs_err_g_kt[l] - pow(obs_mean_g_kt[l], 2), 0.5) / pow(nbin - 1, 0.5);
    }
    obs_mean_rho_s /= nbin;
    obs_err_rho_s /= nbin;
    obs_err_rho_s = pow(obs_err_rho_s - pow(obs_mean_rho_s, 2), 0.5) / pow(nbin - 1, 0.5);
}
