#include "DynamicMeasure.h"
#include "Hubbard.h"

Measure::DynamicMeasure::DynamicMeasure(const int &nbin) {
    this->nbin = nbin;
}

void Measure::DynamicMeasure::resize(const int &_nbin) {
    this->nbin = _nbin;
}

void Measure::DynamicMeasure::initial(const Model::Hubbard &hubbard) {
    this->sign.set_size_of_bin(this->nbin);
    this->superfluid_density.set_size_of_bin(this->nbin);

    this->matsubara_greens.reserve(hubbard.lt);
    this->density_of_states.reserve(hubbard.lt);
    for (int l = 0; l < hubbard.lt; ++l) {
        this->matsubara_greens.emplace_back(this->nbin);
        this->density_of_states.emplace_back(this->nbin);
    }

//    this->current_current_corr.resize(hubbard.ls, hubbard.ls);
//    for (int i = 0; i < hubbard.ls; ++i) {
//        for (int j = 0; j < hubbard.ls; ++j) {
//            this->current_current_corr(i, j).set_size_of_bin(this->nbin);
//        }
//    }
}

void Measure::DynamicMeasure::clear_temporary(const Model::Hubbard &hubbard) {
    this->sign.clear_temporary();
    this->superfluid_density.clear_temporary();
    for (auto &green : this->matsubara_greens) {
        green.clear_temporary();
    }
    for (auto &dos : this->density_of_states) {
        dos.clear_temporary();
    }

//    for (int i = 0; i < hubbard.ls; ++i) {
//        for (int j = 0; j < hubbard.ls; ++j) {
//            this->current_current_corr(i, j).clear_temporary();
//        }
//    }
}

void Measure::DynamicMeasure::time_displaced_measure(const Model::Hubbard &hubbard) {
    this->sign.tmp_value() += hubbard.config_sign;
    ++this->sign;
    for (int l = 0; l < hubbard.lt; ++l) {
        this->measure_matsubara_greens(l, hubbard);
        this->measure_density_of_states(l, hubbard);
    }
    this->measure_superfluid_density(hubbard);
}

void Measure::DynamicMeasure::normalize_stats(const Model::Hubbard &hubbard) {
    this->sign.tmp_value() /= this->sign.counts();
    this->superfluid_density.tmp_value() /= this->superfluid_density.counts() * this->sign.tmp_value();
    for (auto &green : this->matsubara_greens) {
        green.tmp_value() /= green.counts() * this->sign.tmp_value();
    }
    for (auto &dos : this->density_of_states) {
        dos.tmp_value() /= dos.counts() * this->sign.tmp_value();
    }

//    for (int i = 0; i < hubbard.ls; ++i) {
//        for (int j = 0; j < hubbard.ls; ++j) {
//            this->current_current_corr(i, j).tmp_value() /= current_current_corr(i, j).counts() * this->sign.tmp_value();
//        }
//    }
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

//    for (int i = 0; i < hubbard.ls; ++i) {
//        for (int j = 0; j < hubbard.ls; ++j) {
//            this->current_current_corr(i, j).bin_data()[bin] = current_current_corr(j, i).tmp_value();
//        }
//    }
}

void Measure::DynamicMeasure::measure_matsubara_greens(const int &t, const Model::Hubbard &hubbard) {
    assert( t >= 0 && t < hubbard.lt );
    // factor 1/2 comes from two degenerate spin states
    const Eigen::MatrixXd gt0 = ( t == 0 )?
              0.5 * (hubbard.vec_green_tt_up[hubbard.lt-1] + hubbard.vec_green_tt_up[hubbard.lt-1])
            : 0.5 * (hubbard.vec_green_t0_up[t-1] + hubbard.vec_green_t0_up[t-1]);

    // base point i
    for (int xi = 0; xi < hubbard.ll; ++xi) {
        for (int yi = 0; yi < hubbard.ll; ++yi) {
            const int i = xi + hubbard.ll * yi;
            // displacement
            for (int dx = 0; dx < hubbard.ll; ++dx) {
                for (int dy = 0; dy < hubbard.ll; ++dy) {
                    const int j = (xi + dx) % hubbard.ll + hubbard.ll * ((yi + dy) % hubbard.ll);
                    const Eigen::VectorXd r = ( Eigen::VectorXd(2) << dx, dy ).finished();
                    this->matsubara_greens[t].tmp_value() += hubbard.config_sign * cos(-r.dot(this->q)) * gt0(j, i) / hubbard.ls;
                }
            }
        }
    }
    ++this->matsubara_greens[t];
}

void Measure::DynamicMeasure::measure_density_of_states(const int &t, const Model::Hubbard &hubbard) {
    assert( t >= 0 && t < hubbard.lt );
    // spin degenerate model
    const Eigen::MatrixXd gt0 = ( t == 0 )?
              0.5 * (hubbard.vec_green_tt_up[hubbard.lt-1] + hubbard.vec_green_tt_up[hubbard.lt-1])
            : 0.5 * (hubbard.vec_green_t0_up[t-1] + hubbard.vec_green_t0_up[t-1]);
    this->density_of_states[t].tmp_value() += hubbard.config_sign * gt0.trace() / hubbard.ls;
    ++this->density_of_states[t];
}

void Measure::DynamicMeasure::measure_superfluid_density(const Model::Hubbard &hubbard) {
    // momentum qx and qy
    const Eigen::VectorXd qx = ( Eigen::VectorXd(2) << 2 * M_PI / hubbard.ll, 0.0 ).finished();
    const Eigen::VectorXd qy = ( Eigen::VectorXd(2) << 0.0, 2 * M_PI / hubbard.ll ).finished();

    // fourier transformation in time-energy space
    double tmp_rho_s = 0.0;
    Eigen::MatrixXd gt0_up, g0t_up, gtt_up, g00_up;
    Eigen::MatrixXd gt0_dn, g0t_dn, gtt_dn, g00_dn;

    for (int l = 0; l < hubbard.lt; ++l) {
        if ( l == 0 ) {
            gt0_up = hubbard.vec_green_tt_up[hubbard.lt - 1];
            gt0_dn = hubbard.vec_green_tt_dn[hubbard.lt - 1];
            g0t_up = gt0_up;
            gtt_up = gt0_up;
            g00_up = gt0_up;
            g0t_dn = gt0_dn;
            gtt_dn = gt0_dn;
            g00_dn = gt0_dn;
        }
        else {
            gt0_up = hubbard.vec_green_t0_up[l-1];
            g0t_up = hubbard.vec_green_0t_up[l-1];
            gtt_up = hubbard.vec_green_tt_up[l-1];
            g00_up = hubbard.vec_green_tt_up[hubbard.lt - 1];
            gt0_dn = hubbard.vec_green_t0_dn[l-1];
            g0t_dn = hubbard.vec_green_0t_dn[l-1];
            gtt_dn = hubbard.vec_green_tt_dn[l-1];
            g00_dn = hubbard.vec_green_tt_dn[hubbard.lt - 1];
        }

        // space point i is chosen as our base point, which is to be averaged
        for (int xi = 0; xi < hubbard.ll; ++xi) {
            for (int yi = 0; yi < hubbard.ll; ++yi) {
                const int i = xi + hubbard.ll * yi;
                const int ipx = (xi + 1) % hubbard.ll + hubbard.ll * yi;

                // displacement
                for (int dx = 0; dx < hubbard.ll; ++dx) {
                    for (int dy = 0; dy < hubbard.ll; ++dy) {
                        /* for a given site l and time-slice tau
                         * the current-current correlation Jx-Jx: \Gamma_xx (l, \tau) = < jx(l, \tau) * jx(0, 0) > */
                        const int j = (xi + dx) % hubbard.ll + hubbard.ll * ((yi + dy) % hubbard.ll);
                        const int jpx = (xi + dx + 1) % hubbard.ll + hubbard.ll * ((yi + dy) % hubbard.ll);
                        const Eigen::VectorXd r = ( Eigen::VectorXd(2) << dx, dy ).finished();
                        const double factor = hubbard.config_sign * (cos(r.dot(qx)) - cos(r.dot(qy)));

                        tmp_rho_s += hubbard.t * hubbard.t * factor * (
                                - ( gtt_up(j, jpx) - gtt_up(jpx, j) + gtt_dn(j, jpx) - gtt_dn(jpx, j) ) *
                                  ( g00_up(i, ipx) - g00_up(ipx, i) + g00_dn(i, ipx) - g00_dn(ipx, i) )

                                - g0t_up(ipx, jpx) * gt0_up(j, i) - g0t_dn(ipx, jpx) * gt0_dn(j, i)
                                + g0t_up(i, jpx) * gt0_up(j, ipx) + g0t_dn(i, jpx) * gt0_dn(j, ipx)
                                + g0t_up(ipx, j) * gt0_up(jpx, i) + g0t_dn(ipx, j) * gt0_dn(jpx, i)
                                - g0t_up(i, j) * gt0_up(jpx, ipx) - g0t_dn(i, j) * gt0_dn(jpx, ipx)
                        );
                    }
                }
            }
        }
    }
    // average over base point i
    this->superfluid_density.tmp_value() += 0.25 * tmp_rho_s / hubbard.ls / hubbard.ls;
    ++this->superfluid_density;
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

//    for (int i = 0; i < hubbard.ls; ++i) {
//        for (int j = 0; j < hubbard.ls; ++j) {
//            this->current_current_corr(i, j).analyse();
//        }
//    }
}
