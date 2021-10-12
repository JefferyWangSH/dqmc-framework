#include "DynamicMeasure.h"
#include "Hubbard.h"

Measure::DynamicMeasure::DynamicMeasure(const int &nbin) {
    this->nbin = nbin;
}

void Measure::DynamicMeasure::resize(const int &_nbin) {
    this->nbin = _nbin;
}

void Measure::DynamicMeasure::initial(const Model::Hubbard &hubbard) {
    // clear data of previous simulation
    this->clear();

    this->bin_g_kt.reserve(nbin);
    this->bin_rho_s.reserve(nbin);
    this->bin_sign.reserve(nbin);

    this->bin_gt0_up.reserve(nbin);
    this->bin_g0t_up.reserve(nbin);
    this->bin_gtt_up.reserve(nbin);
    this->bin_gt0_dn.reserve(nbin);
    this->bin_g0t_dn.reserve(nbin);
    this->bin_gtt_dn.reserve(nbin);

    for (int bin = 0; bin < this->nbin; ++bin) {
        this->bin_g_kt.emplace_back(hubbard.lt, 0.0);
        this->bin_rho_s.emplace_back(0.0);
        this->bin_sign.emplace_back(0.0);

        this->bin_gt0_up.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        this->bin_g0t_up.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        this->bin_gtt_up.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        this->bin_gt0_dn.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        this->bin_g0t_dn.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        this->bin_gtt_dn.emplace_back(hubbard.lt, Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
    }

    this->mean_rho_s = 0.0;
    this->err_rho_s = 0.0;

    this->mean_sign = 0.0;
    this->err_sign = 0.0;
    this->tmp_sign = 0.0;

    this->mean_g_kt.reserve(hubbard.lt);
    this->err_g_kt.reserve(hubbard.lt);

    this->tmp_gt0_tau_up.reserve(hubbard.lt);
    this->tmp_g0t_tau_up.reserve(hubbard.lt);
    this->tmp_gtt_tau_up.reserve(hubbard.lt);
    this->tmp_gt0_tau_dn.reserve(hubbard.lt);
    this->tmp_g0t_tau_dn.reserve(hubbard.lt);
    this->tmp_gtt_tau_dn.reserve(hubbard.lt);

    for (int l = 0; l < hubbard.lt; ++l) {
        this->mean_g_kt.emplace_back(0.0);
        this->err_g_kt.emplace_back(0.0);

        this->tmp_gt0_tau_up.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        this->tmp_g0t_tau_up.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        this->tmp_gtt_tau_up.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        this->tmp_gt0_tau_dn.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        this->tmp_g0t_tau_dn.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
        this->tmp_gtt_tau_dn.emplace_back(Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls));
    }

    // free redundant memory
    this->bin_g_kt.shrink_to_fit();
    this->mean_g_kt.shrink_to_fit();
    this->err_g_kt.shrink_to_fit();
    this->bin_rho_s.shrink_to_fit();
    this->bin_sign.shrink_to_fit();

    this->bin_gt0_up.shrink_to_fit();
    this->bin_g0t_up.shrink_to_fit();
    this->bin_gtt_up.shrink_to_fit();
    this->bin_gt0_dn.shrink_to_fit();
    this->bin_g0t_dn.shrink_to_fit();
    this->bin_gtt_dn.shrink_to_fit();

    this->tmp_gt0_tau_up.shrink_to_fit();
    this->tmp_g0t_tau_up.shrink_to_fit();
    this->tmp_gtt_tau_up.shrink_to_fit();
    this->tmp_gt0_tau_dn.shrink_to_fit();
    this->tmp_g0t_tau_dn.shrink_to_fit();
    this->tmp_gtt_tau_dn.shrink_to_fit();

    assert( this->bin_gt0_up.size() == this->nbin );
    assert( this->bin_g0t_up.size() == this->nbin );
    assert( this->bin_gtt_up.size() == this->nbin );
    assert( this->bin_gt0_dn.size() == this->nbin );
    assert( this->bin_g0t_dn.size() == this->nbin );
    assert( this->bin_gtt_dn.size() == this->nbin );

    for (int bin = 0; bin < this->nbin; ++bin) {
        this->bin_gt0_up.shrink_to_fit();
        this->bin_g0t_up.shrink_to_fit();
        this->bin_gtt_up.shrink_to_fit();
        this->bin_gt0_dn.shrink_to_fit();
        this->bin_g0t_dn.shrink_to_fit();
        this->bin_gtt_dn.shrink_to_fit();
    }
}

void Measure::DynamicMeasure::clear() {
    this->bin_g_kt.clear();
    this->bin_rho_s.clear();
    this->bin_sign.clear();

    this->bin_gt0_up.clear();
    this->bin_g0t_up.clear();
    this->bin_gtt_up.clear();
    this->bin_gt0_dn.clear();
    this->bin_g0t_dn.clear();
    this->bin_gtt_dn.clear();

    this->mean_g_kt.clear();
    this->err_g_kt.clear();

    this->mean_rho_s = 0.0;
    this->err_rho_s = 0.0;
    this->mean_sign = 0.0;
    this->err_sign = 0.0;

    this->tmp_gt0_tau_up.clear();
    this->tmp_g0t_tau_up.clear();
    this->tmp_gtt_tau_up.clear();
    this->tmp_gt0_tau_dn.clear();
    this->tmp_g0t_tau_dn.clear();
    this->tmp_gtt_tau_dn.clear();
}

void Measure::DynamicMeasure::clear_temporary(const Model::Hubbard &hubbard) {
    assert( this->tmp_gt0_tau_up.size() == hubbard.lt );
    assert( this->tmp_g0t_tau_up.size() == hubbard.lt );
    assert( this->tmp_gtt_tau_up.size() == hubbard.lt );
    assert( this->tmp_gt0_tau_dn.size() == hubbard.lt );
    assert( this->tmp_g0t_tau_dn.size() == hubbard.lt );
    assert( this->tmp_gtt_tau_dn.size() == hubbard.lt );

    this->n_time_displaced = 0;
    this->tmp_sign = 0.0;
    for (int l = 0; l < hubbard.lt; ++l) {
        this->tmp_gt0_tau_up[l].setZero(hubbard.ls, hubbard.ls);
        this->tmp_g0t_tau_up[l].setZero(hubbard.ls, hubbard.ls);
        this->tmp_gtt_tau_up[l].setZero(hubbard.ls, hubbard.ls);
        this->tmp_gt0_tau_dn[l].setZero(hubbard.ls, hubbard.ls);
        this->tmp_g0t_tau_dn[l].setZero(hubbard.ls, hubbard.ls);
        this->tmp_gtt_tau_dn[l].setZero(hubbard.ls, hubbard.ls);
    }
}

void Measure::DynamicMeasure::measure_time_displaced(const Model::Hubbard &hubbard) {
    this->n_time_displaced++;
    for (int l = 0; l < hubbard.lt; ++l) {
        this->tmp_gt0_tau_up[l] += hubbard.config_sign * hubbard.vec_green_t0_up[l];
        this->tmp_g0t_tau_up[l] += hubbard.config_sign * hubbard.vec_green_0t_up[l];
        this->tmp_gtt_tau_up[l] += hubbard.config_sign * hubbard.vec_green_tt_up[l];
        this->tmp_gt0_tau_dn[l] += hubbard.config_sign * hubbard.vec_green_t0_dn[l];
        this->tmp_g0t_tau_dn[l] += hubbard.config_sign * hubbard.vec_green_0t_dn[l];
        this->tmp_gtt_tau_dn[l] += hubbard.config_sign * hubbard.vec_green_tt_dn[l];
    }
    this->tmp_sign += hubbard.config_sign;
}

void Measure::DynamicMeasure::normalizeStats(const Model::Hubbard &hubbard) {
    this->tmp_sign /= this->n_time_displaced;
    for (int l = 0; l < hubbard.lt; ++l) {
        this->tmp_gt0_tau_up[l] /= this->n_time_displaced * this->tmp_sign;
        this->tmp_g0t_tau_up[l] /= this->n_time_displaced * this->tmp_sign;
        this->tmp_gtt_tau_up[l] /= this->n_time_displaced * this->tmp_sign;
        this->tmp_gt0_tau_dn[l] /= this->n_time_displaced * this->tmp_sign;
        this->tmp_g0t_tau_dn[l] /= this->n_time_displaced * this->tmp_sign;
        this->tmp_gtt_tau_dn[l] /= this->n_time_displaced * this->tmp_sign;
    }
}

void Measure::DynamicMeasure::write_Stats_to_bins(const int &bin, const Model::Hubbard &hubbard) {
    // TODO: CHECK IT OUT
    for (int l = 0; l < hubbard.lt; ++l) {
        this->bin_gt0_up[bin][l] = this->tmp_gt0_tau_up[l];
        this->bin_g0t_up[bin][l] = this->tmp_g0t_tau_up[l];
        this->bin_gtt_up[bin][l] = this->tmp_gtt_tau_up[l];
        this->bin_gt0_dn[bin][l] = this->tmp_gt0_tau_dn[l];
        this->bin_g0t_dn[bin][l] = this->tmp_g0t_tau_dn[l];
        this->bin_gtt_dn[bin][l] = this->tmp_gtt_tau_dn[l];
    }
    this->bin_sign[bin] = this->tmp_sign;
}

void Measure::DynamicMeasure::analyse_Dynamical_Corr(const int &bin, const Model::Hubbard &hubbard) {
    for (int l = 0; l < hubbard.lt; ++l) {
        /* factor 0.5 comes from two spin states */
        const Eigen::MatrixXd gt0 = 0.5 * (this->bin_gt0_up[bin][l] + this->bin_gt0_dn[bin][l]);
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
        this->bin_g_kt[bin][l] = tmpFourier;
    }
}

void Measure::DynamicMeasure::analyse_Rho_S(const int &bin, const Model::Hubbard &hubbard) {

    // momentum qx and qy
    const Eigen::VectorXd qx = (Eigen::VectorXd(2) << 2 * M_PI / hubbard.ll, 0).finished();
    const Eigen::VectorXd qy = (Eigen::VectorXd(2) << 0, 2 * M_PI / hubbard.ll).finished();

    double tmp_rho_s = 0.0;

    // fourier transformation in time-energy space
    for (int l = 0; l < hubbard.lt; ++l) {
        double tmp_fourier = 0.0;

        const Eigen::MatrixXd gt0_up = this->bin_gt0_up[bin][l];
        const Eigen::MatrixXd g0t_up = this->bin_g0t_up[bin][l];
        const Eigen::MatrixXd gtt_up = this->bin_gtt_up[bin][l];
        const Eigen::MatrixXd g00_up = this->bin_gtt_up[bin][hubbard.lt - 1];
        const Eigen::MatrixXd gt0_dn = this->bin_gt0_dn[bin][l];
        const Eigen::MatrixXd g0t_dn = this->bin_g0t_dn[bin][l];
        const Eigen::MatrixXd gtt_dn = this->bin_gtt_dn[bin][l];
        const Eigen::MatrixXd g00_dn = this->bin_gtt_dn[bin][hubbard.lt - 1];

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

    this->bin_rho_s[bin] = 0.25 * tmp_rho_s;
}

void Measure::DynamicMeasure::analyse_timeDisplaced_Stats(const Model::Hubbard &hubbard) {
    // analyse time-displaced statistics and record data.

    // clear previous statistics
    for (int l = 0; l < hubbard.lt; ++l) {
        this->mean_g_kt[l] = 0.0;
        this->err_g_kt[l] = 0.0;
    }
    this->mean_rho_s = 0.0;
    this->err_rho_s = 0.0;
    this->mean_sign = 0.0;
    this->err_sign = 0.0;

    // analyse statistics
    for (int bin = 0; bin < this->nbin; ++bin) {
        this->analyse_Dynamical_Corr(bin, hubbard);
        this->analyse_Rho_S(bin, hubbard);
    }

    // sum over bins and calculate mean and error
    for (int bin = 0; bin < this->nbin; ++bin) {
        for (int l = 0; l < hubbard.lt; ++l) {
            this->mean_g_kt[l] += this->bin_g_kt[bin][l];
            this->err_g_kt[l] += this->bin_g_kt[bin][l] * this->bin_g_kt[bin][l];
        }
        this->mean_rho_s += this->bin_rho_s[bin];
        this->err_rho_s += this->bin_rho_s[bin] * this->bin_rho_s[bin];
        this->mean_sign += this->bin_sign[bin];
        this->err_sign += this->bin_sign[bin] * this->bin_sign[bin];
    }

    this->mean_sign /= this->nbin;
    this->err_sign /= this->nbin;
    this->err_sign = pow(this->err_sign - pow(this->mean_sign, 2), 0.5) / pow(this->nbin - 1, 0.5);
    for (int l = 0; l < hubbard.lt; ++l) {
        this->mean_g_kt[l] /= this->nbin;
        this->err_g_kt[l] /= this->nbin;
        this->err_g_kt[l] = pow(this->err_g_kt[l] - pow(this->mean_g_kt[l], 2), 0.5) / pow(this->nbin - 1, 0.5);
    }
    this->mean_rho_s /= this->nbin;
    this->err_rho_s /= this->nbin;
    this->err_rho_s = pow(this->err_rho_s - pow(this->mean_rho_s, 2), 0.5) / pow(this->nbin - 1, 0.5);
}
