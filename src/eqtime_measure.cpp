#include "eqtime_measure.h"
#include "hubbard.h"

Measure::EqtimeMeasure::EqtimeMeasure(const int &nbin) {
    this->nbin = nbin;
}

void Measure::EqtimeMeasure::resize(const int &_nbin) {
    this->nbin = _nbin;
}

void Measure::EqtimeMeasure::initial(const Model::Hubbard &hubbard) {
    this->sign.set_size_of_bin(this->nbin);
    this->filling_number.set_size_of_bin(this->nbin);
    this->double_occupancy.set_size_of_bin(this->nbin);
    this->kinetic_energy.set_size_of_bin(this->nbin);
    this->momentum_distribution.set_size_of_bin(this->nbin);
    this->local_spin_corr.set_size_of_bin(this->nbin);
    this->spin_density_structure_factor.set_size_of_bin(this->nbin);
    this->charge_density_structure_factor.set_size_of_bin(this->nbin);
    this->pairing_corr.reserve(hubbard.ll/2+1);
    for (int i = 0; i < hubbard.ll/2+1; ++i) {
        this->pairing_corr.emplace_back(this->nbin);
    }

    this->sign.set_zero_element(0.0);
    this->filling_number.set_zero_element(0.0);
    this->double_occupancy.set_zero_element(0.0);
    this->kinetic_energy.set_zero_element(0.0);
    this->momentum_distribution.set_zero_element(0.0);
    this->local_spin_corr.set_zero_element(0.0);
    this->spin_density_structure_factor.set_zero_element(0.0);
    this->charge_density_structure_factor.set_zero_element(0.0);
    for (auto &PairingCorr : this->pairing_corr) {
        PairingCorr.set_zero_element(0.0);
    }

    this->sign.allocate();
    this->filling_number.allocate();
    this->double_occupancy.allocate();
    this->kinetic_energy.allocate();
    this->momentum_distribution.allocate();
    this->local_spin_corr.allocate();
    this->spin_density_structure_factor.allocate();
    this->charge_density_structure_factor.allocate();
    for (auto &PairingCorr : this->pairing_corr) {
        PairingCorr.allocate();
    }
}

void Measure::EqtimeMeasure::clear_temporary() {
    this->sign.clear_temporary();
    this->filling_number.clear_temporary();
    this->double_occupancy.clear_temporary();
    this->kinetic_energy.clear_temporary();
    this->momentum_distribution.clear_temporary();
    this->local_spin_corr.clear_temporary();
    this->spin_density_structure_factor.clear_temporary();
    this->charge_density_structure_factor.clear_temporary();
    for (auto &PairingCorr : this->pairing_corr) {
        PairingCorr.clear_temporary();
    }
}

void Measure::EqtimeMeasure::equal_time_measure(const Model::Hubbard &hubbard) {
    this->sign.tmp_value() += hubbard.config_sign;
    ++this->sign;
    for (int t = 0; t < hubbard.lt; ++t) {
        this->measure_filling_number(t, hubbard);
        this->measure_double_occupancy(t, hubbard);
        this->measure_kinetic_energy(t, hubbard);
        this->measure_momentum_distribution(t, hubbard);
        this->measure_local_spin_corr(t, hubbard);
        this->measure_spin_density_structure_factor(t, hubbard);
        this->measure_charge_density_structure_factor(t, hubbard);
        this->measure_pairing_corr(t, hubbard);
    }
}

void Measure::EqtimeMeasure::normalize_stats(const Model::Hubbard &hubbard) {
    // normalized by counting
    this->sign.tmp_value() /= this->sign.counts();
    this->filling_number.tmp_value() /= this->filling_number.counts() * this->sign.tmp_value();
    this->double_occupancy.tmp_value() /= this->double_occupancy.counts() * this->sign.tmp_value();
    this->kinetic_energy.tmp_value() /= this->kinetic_energy.counts() * this->sign.tmp_value();
    this->momentum_distribution.tmp_value() /= this->momentum_distribution.counts() * this->sign.tmp_value();
    this->local_spin_corr.tmp_value() /= this->local_spin_corr.counts() * this->sign.tmp_value();
    this->spin_density_structure_factor.tmp_value() /= this->spin_density_structure_factor.counts() * this->sign.tmp_value();
    this->charge_density_structure_factor.tmp_value() /= this->charge_density_structure_factor.counts() * this->sign.tmp_value();
    for (auto &PairingCorr : this->pairing_corr) {
        PairingCorr.tmp_value() /= PairingCorr.counts() * this->sign.tmp_value();
    }
}

void Measure::EqtimeMeasure::write_stats_to_bins(int bin) {
    this->sign.bin_data()[bin] = this->sign.tmp_value();
    this->filling_number.bin_data()[bin] = this->filling_number.tmp_value();
    this->double_occupancy.bin_data()[bin] = this->double_occupancy.tmp_value();
    this->kinetic_energy.bin_data()[bin] = this->kinetic_energy.tmp_value();
    this->momentum_distribution.bin_data()[bin] = this->momentum_distribution.tmp_value();
    this->local_spin_corr.bin_data()[bin] = this->local_spin_corr.tmp_value();
    this->spin_density_structure_factor.bin_data()[bin] = this->spin_density_structure_factor.tmp_value();
    this->charge_density_structure_factor.bin_data()[bin] = this->charge_density_structure_factor.tmp_value();
    for (auto &PairingCorr : this->pairing_corr) {
        PairingCorr.bin_data()[bin] = PairingCorr.tmp_value();
    }
}

void Measure::EqtimeMeasure::measure_filling_number(const int &t, const Model::Hubbard &hubbard) {
    assert( t >= 0 && t < hubbard.lt );
    const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
    const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];
    this->filling_number.tmp_value() += hubbard.config_sign * (2 - (gu.trace()+gd.trace())/hubbard.ls);
    ++this->filling_number;
}

void Measure::EqtimeMeasure::measure_double_occupancy(const int &t, const Model::Hubbard &hubbard) {
    assert( t >= 0 && t < hubbard.lt );
    const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
    const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

    double tmp_double_occu = 0.0;
    for (int i = 0; i < hubbard.ls; ++i) {
        tmp_double_occu += hubbard.config_sign * (1 - gu(i, i)) * (1 - gd(i, i));
    }
    this->double_occupancy.tmp_value() += tmp_double_occu / hubbard.ls;
    ++this->double_occupancy;
}

void Measure::EqtimeMeasure::measure_kinetic_energy(const int &t, const Model::Hubbard &hubbard) {
    assert( t >= 0 && t < hubbard.lt );
    const int ll = hubbard.ll;
    const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
    const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

    double tmp_kinetic_energy = 0.0;
    for (int x = 0; x < ll; ++x) {
        for (int y = 0; y < ll; ++y) {
            tmp_kinetic_energy += 2 * hubbard.t * hubbard.config_sign *
                ( gu(x + ll*y, ((x+1)%ll) + ll*y) + gu(x + ll*y, x + ll*((y+1)%ll))
                + gd(x + ll*y, ((x+1)%ll) + ll*y) + gd(x + ll*y, x + ll*((y+1)%ll)) );
        }
    }
    this->kinetic_energy.tmp_value() += tmp_kinetic_energy / hubbard.ls;
    ++this->kinetic_energy;
}

void Measure::EqtimeMeasure::measure_momentum_distribution(const int &t, const Model::Hubbard &hubbard) {
    assert( t >= 0 && t < hubbard.lt );
    const int ll = hubbard.ll;
    const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
    const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

    double tmp_momentum_dist = 0.0;
    // base point
    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            const int i = xi + ll * yi;
            // displacement
            for (int dx = 0; dx < ll; ++dx) {
                for (int dy = 0; dy < ll; ++dy) {
                    const int j = (xi + dx) % ll + ll * ((yi + dy) % ll);
                    const Eigen::VectorXd r = ( Eigen::VectorXd(2) << dx, dy ).finished();
                    tmp_momentum_dist += cos(-r.dot(this->q)) * (gu(j, i) + gd(j, i));
                }
            }
        }
    }
    this->momentum_distribution.tmp_value() += hubbard.config_sign * (1 - 0.5 * tmp_momentum_dist / hubbard.ls);
    ++this->momentum_distribution;
}

void Measure::EqtimeMeasure::measure_local_spin_corr(const int &t, const Model::Hubbard &hubbard) {
    assert( t >= 0 && t < hubbard.lt );
    const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
    const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

    double tmp_local_spin_corr = 0.0;
    for (int i = 0; i < hubbard.ls; ++i) {
        tmp_local_spin_corr += hubbard.config_sign * ( gu(i, i) + gd(i, i) - 2 * gu(i, i) * gd(i, i) );
    }
    this->local_spin_corr.tmp_value() += tmp_local_spin_corr / hubbard.ls;
    ++this->local_spin_corr;
}

void Measure::EqtimeMeasure::measure_spin_density_structure_factor(const int &t, const Model::Hubbard &hubbard) {
    assert( t >= 0 && t < hubbard.lt );
    const int ll = hubbard.ll;
    const int ls = hubbard.ls;
    const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
    const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

    // g(i,j)  = < c_i * c^+_j >
    // gc(i,j) = < c^+_i * c_j >
    const Eigen::MatrixXd guc = Eigen::MatrixXd::Identity(ls, ls) - gu.transpose();
    const Eigen::MatrixXd gdc = Eigen::MatrixXd::Identity(ls, ls) - gd.transpose();

    // loop for site i and average
    double tmp_sdw = 0.0;
    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            const int i = xi + ll * yi;
            // displacement
            for (int dx = 0; dx < ll; ++dx) {
                for (int dy = 0; dy < ll; ++dy) {
                    const int j = (xi + dx) % ll + ll * ((yi + dy) % ll);
                    const Eigen::VectorXd r = ( Eigen::VectorXd(2) << dx, dy ).finished();
                    const double factor = hubbard.config_sign * cos(-r.dot(this->q));
                    // factor 1/4 comes from spin 1/2
                    tmp_sdw += 0.25 * factor * ( + guc(i, i) * guc(j, j) + guc(i, j) * gu(i, j)
                                                 + gdc(i, i) * gdc(j, j) + gdc(i, j) * gd(i, j)
                                                 - gdc(i, i) * guc(j, j) - guc(i, i) * gdc(j, j) );
                }
            }
        }
    }
    this->spin_density_structure_factor.tmp_value() += tmp_sdw / hubbard.ls / hubbard.ls;
    ++this->spin_density_structure_factor;
}

void Measure::EqtimeMeasure::measure_charge_density_structure_factor(const int &t, const Model::Hubbard &hubbard) {
    assert( t >= 0 && t < hubbard.lt );
    const int ll = hubbard.ll;
    const int ls = hubbard.ls;
    const Eigen::MatrixXd gu = (*hubbard.vec_green_tt_up)[t];
    const Eigen::MatrixXd gd = (*hubbard.vec_green_tt_dn)[t];

    // g(i,j)  = < c_i * c^+_j >
    // gc(i,j) = < c^+_i * c_j >
    const Eigen::MatrixXd guc = Eigen::MatrixXd::Identity(ls, ls) - gu.transpose();
    const Eigen::MatrixXd gdc = Eigen::MatrixXd::Identity(ls, ls) - gd.transpose();

    // loop for site i and average
    double tmp_cdw = 0.0;
    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            const int i = xi + ll * yi;
            // displacement
            for (int dx = 0; dx < ll; ++dx) {
                for (int dy = 0; dy < ll; ++dy) {
                    const int j = (xi + dx) % ll + ll * ((yi + dy) % ll);
                    const Eigen::VectorXd r = ( Eigen::VectorXd(2) << dx, dy ).finished();
                    const double factor = hubbard.config_sign * cos(-r.dot(this->q));
                    tmp_cdw += factor * ( + guc(i, i) * guc(j, j) + guc(i, j) * gu(i, j)
                                          + gdc(i, i) * gdc(j, j) + gdc(i, j) * gd(i, j)
                                          + gdc(i, i) * guc(j, j) + guc(i, i) * gdc(j, j) );
                }
            }
        }
    }
    this->charge_density_structure_factor.tmp_value() += tmp_cdw / hubbard.ls / hubbard.ls;
    ++this->charge_density_structure_factor;

}

void Measure::EqtimeMeasure::measure_pairing_corr(const int &t, const Model::Hubbard &hubbard) {
    assert( t >= 0 && t < hubbard.lt );
    // The s-wave pairing order parameter is defined as \Delta_{i} = c_up_{i} * c_dn_{i}
    // Accordingly, the space correlation of Cooper pair reads as follows
    //   < \Delta^\dagger_{i} * \Delta_{j} >
    //        = < c_up_{i}^+ * c_up_{j} > * < c_dn_{i}^+ * c_dn_{j} >
    //        = ( \delta_{ij} - G_up(tau, tau)_{ji} ) * ( \delta_{ij} - G_dn(tau, tau)_{ji} )

    // g(i,j)  = < c_i * c^+_j >
    // gc(i,j) = < c^+_i * c_j >
    const int ll = hubbard.ll;
    const int ls = hubbard.ls;
    const Eigen::MatrixXd guc = Eigen::MatrixXd::Identity(ls, ls) - (*hubbard.vec_green_tt_up)[t].transpose();
    const Eigen::MatrixXd gdc = Eigen::MatrixXd::Identity(ls, ls) - (*hubbard.vec_green_tt_dn)[t].transpose();

    // number of independent correlation length = ll/2 + 1 due to the periodical boundary condition
    // correlation length l1, l2 satisfying l1 + l2 = ll are identical
    const int num_of_corr = ll / 2 + 1;

    // select base point i
    for (int xi = 0; xi < ll; ++xi) {
        for (int yi = 0; yi < ll; ++yi) {
            const int i = xi + ll * yi;

            // correlations along x direction, including onsite cases (i=j)
            for (int dx = 0; dx < ll; ++dx) {
                const int j = (xi + dx) % ll + ll * yi;
                const int corr_len = (dx < num_of_corr)? dx : ll-dx;
                this->pairing_corr[corr_len].tmp_value() += hubbard.config_sign * guc(i, j) * gdc(i, j);
                ++this->pairing_corr[corr_len];
            }

            // correlations along y direction
            for (int dy = 1; dy < ll; ++dy) {
                const int j = xi + ll * ((yi + dy) % ll);
                const int corr_len = (dy < num_of_corr)? dy : ll-dy;
                this->pairing_corr[corr_len].tmp_value() += hubbard.config_sign * guc(i, j) * gdc(i, j);
                ++this->pairing_corr[corr_len];
            }
        }
    }
}

void Measure::EqtimeMeasure::analyse_stats(const Model::Hubbard &hubbard) {
    this->sign.analyse();
    this->filling_number.analyse();
    this->double_occupancy.analyse();
    this->kinetic_energy.analyse();
    this->momentum_distribution.analyse();
    this->local_spin_corr.analyse();
    this->spin_density_structure_factor.analyse();
    this->charge_density_structure_factor.analyse();
    for (auto &PairingCorr : this->pairing_corr) {
        PairingCorr.analyse();
    }
}
