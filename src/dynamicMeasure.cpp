#include "dynamicMeasure.h"
#include "hubbard.h"

void measure::dynamicMeasure::resize(const int &nbin) {
    this->nbin = nbin;
}

void measure::dynamicMeasure::initial(const Hubbard &hubbard) {
    obs_bin_g_kt.reserve(nbin);
    obs_bin_rho_s.reserve(nbin);

    for (int bin = 0; bin < nbin; ++bin) {
        obs_bin_g_kt.emplace_back(hubbard.lt, 0.0);
        obs_bin_rho_s.emplace_back(0.0);
        for (int l = 0; l < hubbard.lt; ++l) {
            obs_bin_gt0_up[bin][l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
            obs_bin_g0t_up[bin][l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
            obs_bin_gtt_up[bin][l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
            obs_bin_gt0_dn[bin][l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
            obs_bin_g0t_dn[bin][l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
            obs_bin_gtt_dn[bin][l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
        }
    }

    obs_mean_rho_s = 0.0;
    obs_err_rho_s = 0.0;
    obs_mean_g_kt.reserve(hubbard.lt);
    obs_err_g_kt.reserve(hubbard.lt);
    for (int l = 0; l < hubbard.lt; ++l) {
        obs_mean_g_kt.emplace_back(0.0);
        obs_err_g_kt.emplace_back(0.0);
        tmp_gt0_tau_up[l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
        tmp_g0t_tau_up[l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
        tmp_gtt_tau_up[l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
        tmp_gt0_tau_dn[l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
        tmp_g0t_tau_dn[l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
        tmp_gtt_tau_dn[l] = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
    }
}

void measure::dynamicMeasure::clear(const Hubbard &hubbard) {
    n_time_displaced = 0;
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
        tmp_gt0_tau_up[l] += hubbard.vecGreen_t0_up[l];
        tmp_g0t_tau_up[l] += hubbard.vecGreen_0t_up[l];
        tmp_gtt_tau_up[l] += hubbard.vecGreenU[l];
        tmp_gt0_tau_dn[l] += hubbard.vecGreen_t0_dn[l];
        tmp_g0t_tau_dn[l] += hubbard.vecGreen_0t_dn[l];
        tmp_gtt_tau_dn[l] += hubbard.vecGreenD[l];
    }
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
}

void measure::dynamicMeasure::write_Stats_to_bins(const int &bin, const Hubbard &hubbard) {
    for (int l = 0; l < hubbard.lt; ++l) {
        obs_bin_gt0_up[bin][l] = tmp_gt0_tau_up[l];
        obs_bin_g0t_up[bin][l] = tmp_g0t_tau_up[l];
        obs_bin_gtt_up[bin][l] = tmp_gtt_tau_up[l];
        obs_bin_gt0_dn[bin][l] = tmp_gt0_tau_dn[l];
        obs_bin_g0t_dn[bin][l] = tmp_g0t_tau_dn[l];
        obs_bin_gtt_dn[bin][l] = tmp_gtt_tau_dn[l];
    }
}

void measure::dynamicMeasure::analyse_Dynamical_Corr(const int &bin, const Hubbard &hubbard) {
    for (int l = 0; l < hubbard.lt; ++l) {
        /* factor 2 comes from two spin states */
        const Eigen::MatrixXd gt0 = (obs_bin_gt0_up[bin][l] + obs_bin_gt0_dn[bin][l]) / 2;
        double tmpFourier = 0.0;

        for (int xi = 0; xi < hubbard.ll; ++xi) {
            for (int yi = 0; yi < hubbard.ll; ++yi) {
                for (int xj = 0; xj < hubbard.ll; ++xj) {
                    for (int yj = 0; yj < hubbard.ll; ++yj) {
                        const int i = xi + hubbard.ll * yi;
                        const int j = xj + hubbard.ll * yj;
                        const Eigen::VectorXd r = (Eigen::VectorXd(2) << (xi-xj), (yi-yj)).finished();
                        // TODO: check complex or real and check factor
                        tmpFourier += cos(-r.dot(q)) * gt0(j, i) / hubbard.ls;
                    }
                }
            }
        }
        obs_bin_g_kt[bin][l] = tmpFourier;
    }
}

void measure::dynamicMeasure::analyse_Rho_S(const int &bin, const Hubbard &hubbard) {

    // datum point in space
    const int i0 = 0 + hubbard.ll * 0;
    const int j0 = 1 + hubbard.ll * 0;

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

        for (int x = 0; x < hubbard.ll; ++x) {
            for (int y = 0; y < hubbard.ll; ++y) {
                /* for a given site l(x, y) and time-slice tau
                 * the current-current correlation: \Gamma_xx (l, \tau) = < jx(l, \tau) * jx(0, 0) > */
                const int i = x + hubbard.ll * y;
                const int j = (x + 1) % hubbard.ll + hubbard.ll * y;
                const Eigen::VectorXd r = (Eigen::VectorXd(2) << x, y).finished();
                const double factor = cos(r.dot(qx)) - cos(r.dot(qy));

                tmp_fourier += - hubbard.t * hubbard.t * factor * (
                            + (gtt_up(i, j) - gtt_up(j, i)) * (g00_up(i0, j0) - g00_up(j0, i0))
                            + (gtt_up(i, j) - gtt_up(j, i)) * (g00_dn(i0, j0) - g00_up(j0, i0))
                            + (gtt_dn(i, j) - gtt_up(j, i)) * (g00_up(i0, j0) - g00_up(j0, i0))
                            + (gtt_dn(i, j) - gtt_up(j, i)) * (g00_dn(i0, j0) - g00_up(j0, i0))
                            - g0t_up(i0, j) * gt0_up(i, j0) - g0t_up(j0, i) * gt0_up(j, i0)
                            + g0t_up(j0, j) * gt0_up(i, i0) + g0t_up(i0, i) * gt0_up(j ,j0)
                            - g0t_dn(i0, j) * gt0_dn(i, j0) - g0t_dn(j0, i) * gt0_dn(j, i0)
                            + g0t_dn(j0, j) * gt0_dn(i, i0) + g0t_dn(i0, i) * gt0_dn(j ,j0)
                            );
            }
        }
        tmp_fourier /= hubbard.ls;
        tmp_rho_s += tmp_fourier;       // todo: check factor beta / lt
    }

    obs_bin_rho_s[bin] = tmp_rho_s / 4;
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
    }

    for (int l = 0; l < hubbard.lt; ++l) {
        obs_mean_g_kt[l] /= nbin;
        obs_err_g_kt[l] /= nbin;
        obs_err_g_kt[l] = pow(obs_err_g_kt[l] - pow(obs_mean_g_kt[l], 2), 0.5) / pow(nbin - 1, 0.5);
    }
    obs_mean_rho_s /= nbin;
    obs_err_rho_s /= nbin;
    obs_err_rho_s = pow(obs_err_rho_s - pow(obs_mean_rho_s, 2), 0.5) / pow(nbin - 1, 0.5);
}
