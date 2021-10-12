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

    this->bin_gtt_up.reserve(this->nbin);
    this->bin_gtt_dn.reserve(this->nbin);
    for (int bin = 0; bin < this->nbin; ++bin) {
        this->bin_gtt_up.emplace_back(hubbard.ls, hubbard.ls);
        this->bin_gtt_dn.emplace_back(hubbard.ls, hubbard.ls);
    }

    this->tmp_gtt_up = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
    this->tmp_gtt_dn = Eigen::MatrixXd::Zero(hubbard.ls, hubbard.ls);
}

void Measure::EqtimeMeasure::clear_temporary() {
    this->n_equal_time = 0;
    this->tmp_sign = 0.0;
    this->tmp_gtt_up.setZero();
    this->tmp_gtt_dn.setZero();
}

void Measure::EqtimeMeasure::measure_equal_time_greens(const Model::Hubbard &hubbard) {
    this->tmp_sign += hubbard.config_sign;
    this->n_equal_time++;

    for (int l = 0; l < hubbard.lt; ++l) {
        this->tmp_gtt_up += hubbard.config_sign * hubbard.vec_green_tt_up[l];
        this->tmp_gtt_dn += hubbard.config_sign * hubbard.vec_green_tt_dn[l];
    }
}

void Measure::EqtimeMeasure::normalizeStats(const Model::Hubbard &hubbard) {
    this->tmp_sign /= this->n_equal_time;
    this->tmp_gtt_up /= this->n_equal_time * hubbard.lt * this->tmp_sign;
    this->tmp_gtt_dn /= this->n_equal_time * hubbard.lt * this->tmp_sign;
}

void Measure::EqtimeMeasure::write_Stats_to_bins(int bin) {
    this->sign.bin_data()[bin] = this->tmp_sign;
    this->bin_gtt_up[bin] = this->tmp_gtt_up;
    this->bin_gtt_dn[bin] = this->tmp_gtt_dn;
}

void Measure::EqtimeMeasure::analyse_double_occu(const int &bin, const Model::Hubbard &hubbard) {
    const Eigen::MatrixXd gu = this->bin_gtt_up[bin];
    const Eigen::MatrixXd gd = this->bin_gtt_dn[bin];

    double tmp_double_occu = 0.0;
    for (int i = 0; i < hubbard.ls; ++i) {
        tmp_double_occu += (1 - gu(i, i)) * (1 - gd(i, i));
    }
    this->double_occu.bin_data()[bin] = tmp_double_occu / hubbard.ls;
}

void Measure::EqtimeMeasure::analyse_kinetic_energy(const int &bin, const Model::Hubbard &hubbard) {
    const int ll = hubbard.ll;
    const Eigen::MatrixXd gu = this->bin_gtt_up[bin];
    const Eigen::MatrixXd gd = this->bin_gtt_dn[bin];

    double tmp_kinetic_energy = 0.0;
    for (int x = 0; x < ll; ++x) {
        for (int y = 0; y < ll; ++y) {
            tmp_kinetic_energy += 2 * hubbard.t * (gu(x + ll*y, ((x+1)%ll) + ll*y) + gu(x + ll*y, x + ll*((y+1)%ll)))
                                + 2 * hubbard.t * (gd(x + ll*y, ((x+1)%ll) + ll*y) + gd(x + ll*y, x + ll*((y+1)%ll)));
        }
    }
    this->kinetic_energy.bin_data()[bin] = tmp_kinetic_energy / hubbard.ls;
}

void Measure::EqtimeMeasure::analyse_electron_density(const int &bin, const Model::Hubbard &hubbard) {
    const int ll = hubbard.ll;
    const Eigen::MatrixXd gu = this->bin_gtt_up[bin];
    const Eigen::MatrixXd gd = this->bin_gtt_dn[bin];

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
    tmp_electron_density = (1 - 0.5 * tmp_electron_density / hubbard.ls);
    this->electron_density.bin_data()[bin] = tmp_electron_density;
}

void Measure::EqtimeMeasure::analyse_local_corr(const int &bin, const Model::Hubbard &hubbard) {
    const Eigen::MatrixXd gu = this->bin_gtt_up[bin];
    const Eigen::MatrixXd gd = this->bin_gtt_dn[bin];

    double tmp_local_corr = 0.0;
    for (int i = 0; i < hubbard.ls; ++i) {
        tmp_local_corr += gu(i, i) + gd(i, i) - 2 * gu(i, i) * gd(i, i);
    }
    this->local_corr.bin_data()[bin] = tmp_local_corr / hubbard.ls;
}

void Measure::EqtimeMeasure::analyse_AFM_factor(const int &bin, const Model::Hubbard &hubbard) {
    const int ll = hubbard.ll;
    const int ls = hubbard.ls;
    const Eigen::MatrixXd gu = this->bin_gtt_up[bin];
    const Eigen::MatrixXd gd = this->bin_gtt_dn[bin];

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
    double tmp_AFM_factor = 0.0;
    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            for (int xj = 0; xj < ll; ++xj) {
                for (int yj = 0; yj < ll; ++yj) {
                    const int i = xi + ll * yi;
                    const int j = xj + ll * yj;
                    const Eigen::VectorXd r = (Eigen::VectorXd(2) << (xi - xj), (yi - yj)).finished();
                    const double factor = cos(-r.dot(this->q));
                    // factor 1/4 comes from spin 1/2
                    tmp_AFM_factor += 0.25 * factor * ( + guc(i, i) * guc(j, j) + guc(i, j) * gu(i, j)
                                                        + gdc(i, i) * gdc(j, j) + gdc(i, j) * gd(i, j)
                                                        - gdc(i, i) * guc(j, j) - guc(i, i) * gdc(j, j) );
                }
            }
        }
    }
    this->AFM_factor.bin_data()[bin] = tmp_AFM_factor / (ls * ls);
}

void Measure::EqtimeMeasure::analyse_Cooper_corr(const int &bin, const Model::Hubbard &hubbard) {

}

void Measure::EqtimeMeasure::analyseStats(const Model::Hubbard &hubbard) {
    for (int bin = 0; bin < this->nbin; ++bin) {
        this->analyse_double_occu(bin, hubbard);
        this->analyse_kinetic_energy(bin, hubbard);
        this->analyse_electron_density(bin, hubbard);
        this->analyse_local_corr(bin, hubbard);
        this->analyse_AFM_factor(bin, hubbard);
        this->analyse_Cooper_corr(bin, hubbard);
    }

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
