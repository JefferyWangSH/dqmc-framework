#include "EqtimeMeasure.h"
#include "MeasureData.h"
#include "Hubbard.h"

Measure::EqtimeMeasure::EqtimeMeasure(const int &nbin) {
    this->nbin = nbin;
}

void Measure::EqtimeMeasure::resize(const int &_nbin) {
    this->nbin = _nbin;
}

void Measure::EqtimeMeasure::initial(const Model::Hubbard &hubbard) {
    this->n_equal_time = 0;

    this->double_occu.set_size_of_bin(this->nbin);
    this->kinetic_energy.set_size_of_bin(this->nbin);
    this->electron_density.set_size_of_bin(this->nbin);
    this->local_corr.set_size_of_bin(this->nbin);
    this->AFM_factor.set_size_of_bin(this->nbin);
    this->sign.set_size_of_bin(this->nbin);

    this->cooper_corr.reserve(hubbard.ll);
    for (int l = 0; l < hubbard.ll; ++l) {
        this->cooper_corr.emplace_back(this->nbin);
    }
}

void Measure::EqtimeMeasure::clear_temporary() {
    this->n_equal_time = 0;
    this->sign.clear_temporary();
    this->double_occu.clear_temporary();
    this->kinetic_energy.clear_temporary();
    this->electron_density.clear_temporary();
    this->local_corr.clear_temporary();
    this->AFM_factor.clear_temporary();
    for (auto &CooperCorr : this->cooper_corr) {
        CooperCorr.clear_temporary();
    }
}

void Measure::EqtimeMeasure::equal_time_measure(const Model::Hubbard &hubbard) {
    this->sign.tmp_value() += hubbard.config_sign;
    for (int t = 0; t < hubbard.lt; ++t) {
        this->measure_double_occu(t, hubbard);
        this->measure_kinetic_energy(t, hubbard);
        this->measure_electron_density(t, hubbard);
        this->measure_local_corr(t, hubbard);
        this->measure_AFM_factor(t, hubbard);
        this->measure_Cooper_corr(t, hubbard);
    }
    this->n_equal_time++;
}

void Measure::EqtimeMeasure::normalize_stats(const Model::Hubbard &hubbard) {
    // normalized by counting
    this->sign.tmp_value() /= this->n_equal_time;
    this->double_occu.tmp_value() /= this->n_equal_time * hubbard.lt;
    this->kinetic_energy.tmp_value() /= this->n_equal_time * hubbard.lt;
    this->electron_density.tmp_value() /= this->n_equal_time * hubbard.lt;
    this->local_corr.tmp_value() /= this->n_equal_time * hubbard.lt;
    this->AFM_factor.tmp_value() /= this->n_equal_time * hubbard.lt;
    for (auto &CooperCorr : this->cooper_corr) {
        CooperCorr.tmp_value() /= this->n_equal_time * hubbard.lt;
    }

    // normalized by pre-factor of nature
    this->double_occu.tmp_value() /= hubbard.ls * this->sign.tmp_value();
    this->kinetic_energy.tmp_value() /= hubbard.ls * this->sign.tmp_value();
    this->electron_density.tmp_value() /= this->sign.tmp_value();
    this->local_corr.tmp_value() /= hubbard.ls * this->sign.tmp_value();
    this->AFM_factor.tmp_value() /= hubbard.ls * hubbard.ls * this->sign.tmp_value();
//    for (auto &CooperCorr : this->cooper_corr) {
//        CooperCorr.tmp_value() /= this->n_equal_time * this->sign.tmp_value();
//    }

}

void Measure::EqtimeMeasure::write_stats_to_bins(int bin) {
    this->sign.bin_data()[bin] = this->sign.tmp_value();
    this->double_occu.bin_data()[bin] = this->double_occu.tmp_value();
    this->kinetic_energy.bin_data()[bin] = this->kinetic_energy.tmp_value();
    this->electron_density.bin_data()[bin] = this->electron_density.tmp_value();
    this->local_corr.bin_data()[bin] = this->local_corr.tmp_value();
    this->AFM_factor.bin_data()[bin] = this->AFM_factor.tmp_value();
    for (auto &CooperCorr : this->cooper_corr) {
        CooperCorr.bin_data()[bin] = CooperCorr.tmp_value();
    }
}

void Measure::EqtimeMeasure::measure_double_occu(const int &t, const Model::Hubbard &hubbard) {
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];

    for (int i = 0; i < hubbard.ls; ++i) {
        this->double_occu.tmp_value() += hubbard.config_sign * (1 - gu(i, i)) * (1 - gd(i, i));
    }
}

void Measure::EqtimeMeasure::measure_kinetic_energy(const int &t, const Model::Hubbard &hubbard) {
    const int ll = hubbard.ll;
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];

    for (int x = 0; x < ll; ++x) {
        for (int y = 0; y < ll; ++y) {
            this->kinetic_energy.tmp_value() += 2 * hubbard.t * hubbard.config_sign *
                    ( gu(x + ll*y, ((x+1)%ll) + ll*y) + gu(x + ll*y, x + ll*((y+1)%ll))
                    + gd(x + ll*y, ((x+1)%ll) + ll*y) + gd(x + ll*y, x + ll*((y+1)%ll)) );
        }
    }
}

void Measure::EqtimeMeasure::measure_electron_density(const int &t, const Model::Hubbard &hubbard) {
    const int ll = hubbard.ll;
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];

    double tmp_electron_density = 0.0;
    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            for (int xj = 0; xj < ll; ++xj) {
                for (int yj = 0; yj < ll; ++yj) {
                    const int i = xi + ll * yi;
                    const int j = xj + ll * yj;
                    const Eigen::VectorXd r = ( Eigen::VectorXd(2) << (xi - xj), (yi - yj) ).finished();
                    tmp_electron_density += cos(-r.dot(this->q)) * (gu(j, i) + gd(j, i));
                }
            }
        }
    }
    this->electron_density.tmp_value() += hubbard.config_sign * (1 - 0.5 * tmp_electron_density / hubbard.ls);
}

void Measure::EqtimeMeasure::measure_local_corr(const int &t, const Model::Hubbard &hubbard) {
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];

    for (int i = 0; i < hubbard.ls; ++i) {
        this->local_corr.tmp_value() += hubbard.config_sign * ( gu(i, i) + gd(i, i) - 2 * gu(i, i) * gd(i, i) );
    }
}

void Measure::EqtimeMeasure::measure_AFM_factor(const int &t, const Model::Hubbard &hubbard) {
    const int ll = hubbard.ll;
    const int ls = hubbard.ls;
    const Eigen::MatrixXd gu = hubbard.vec_green_tt_up[t];
    const Eigen::MatrixXd gd = hubbard.vec_green_tt_dn[t];

    // Definition:
    //   gu(i,j)  = < c_i * c^+_j >
    //   guc(i,j) = < c^+_i * c_j >
    Eigen::MatrixXd guc = Eigen::MatrixXd::Identity(ls, ls);
    Eigen::MatrixXd gdc = Eigen::MatrixXd::Identity(ls, ls);
    for (int i = 0; i < ls; ++i) {
        for (int j = 0; j < ls; ++j) {
            guc(j, i) = - gu(i, j);
            gdc(j, i) = - gd(i, j);
        }
        guc(i, i)++;
        gdc(i, i)++;
    }

    // loop for site i and j
    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            for (int xj = 0; xj < ll; ++xj) {
                for (int yj = 0; yj < ll; ++yj) {
                    const int i = xi + ll * yi;
                    const int j = xj + ll * yj;
                    const Eigen::VectorXd r = (Eigen::VectorXd(2) << (xi - xj), (yi - yj)).finished();
                    const double factor = cos(-r.dot(this->q));
                    // factor 1/4 comes from spin 1/2
                    this->AFM_factor.tmp_value() += 0.25 * factor * hubbard.config_sign *
                            ( + guc(i, i) * guc(j, j) + guc(i, j) * gu(i, j)
                              + gdc(i, i) * gdc(j, j) + gdc(i, j) * gd(i, j)
                              - gdc(i, i) * guc(j, j) - guc(i, i) * gdc(j, j) );
                }
            }
        }
    }
}

void Measure::EqtimeMeasure::measure_Cooper_corr(const int &t, const Model::Hubbard &hubbard) {

}

void Measure::EqtimeMeasure::analyse_stats(const Model::Hubbard &hubbard) {
    this->sign.analyse();
    this->double_occu.analyse();
    this->kinetic_energy.analyse();
    this->electron_density.analyse();
    this->local_corr.analyse();
    this->AFM_factor.analyse();
    for (int i = 0; i < hubbard.ll; ++i) {
        this->cooper_corr[i].analyse();
    }
}
