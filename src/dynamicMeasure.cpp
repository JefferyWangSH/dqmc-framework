#include "dynamicMeasure.h"
#include "hubbard.h"

void measure::dynamicMeasure::resize(const int &nbin) {
    this->nbin = nbin;
}

void measure::dynamicMeasure::initial(const Hubbard &hubbard) {
    for (int bin = 0; bin < nbin; ++bin) {
        for (int l = 0; l < hubbard.lt; ++l) {
            obs_bin_gt0[bin][l] = matXd::Zero(hubbard.ls, hubbard.ls);
        }
    }

    for (int l = 0; l < hubbard.lt; ++l) {
        obs_mean_gt0_k[l] = 0.0;
        obs_err_gt0_k[l] = 0.0;
        vec_gt0_tau[l] = matXd::Identity(hubbard.ls, hubbard.ls);
    }
}

void measure::dynamicMeasure::clear(const Hubbard &hubb) {
    n_time_displaced = 0;
    for (int l = 0; l < hubb.lt; ++l) {
        vec_gt0_tau[l] = matXd::Zero(hubb.ls, hubb.ls);
    }
}

void measure::dynamicMeasure::measure_time_displaced(const Hubbard &hubbard) {
    n_time_displaced++;
    for (int l = 0; l < hubbard.lt; ++l) {
        /* factor 2 comes from two spin states */
        vec_gt0_tau[l] += (hubbard.vecGreen_t0_up[l] + hubbard.vecGreen_t0_dn[l]) / 2;
    }
}

void measure::dynamicMeasure::normalizeStats(const Hubbard &hubbard) {
    for (int l = 0; l < hubbard.lt; ++l) {
        vec_gt0_tau[l] /= n_time_displaced;
    }
}

void measure::dynamicMeasure::write_Stats_to_bins(const int &bin, const Hubbard &hubbard) {
    for (int l = 0; l < hubbard.lt; ++l) {
        obs_bin_gt0[bin][l] = vec_gt0_tau[l];
    }
}

void measure::dynamicMeasure::analyse_timeDisplaced_Stats(const Hubbard &hubbard) {
    // calculate the Matsubara green function in momentum space and record data.

    for (int l = 0; l < hubbard.lt; ++l) {
        // clear previous statistics
        obs_mean_gt0_k[l] = 0.0;
        obs_err_gt0_k[l] = 0.0;

        for (int bin = 0; bin < nbin; ++bin) {
            const matXd Matsubara = obs_bin_gt0[bin][l];
            double tmpFourier = 0.0;

            for (int xi = 0; xi < hubbard.ll; ++xi) {
                for (int yi = 0; yi < hubbard.ll; ++yi) {
                    for (int xj = 0; xj < hubbard.ll; ++xj) {
                        for (int yj = 0; yj < hubbard.ll; ++yj) {
                            const int i = xi + hubbard.ll * yi;
                            const int j = xj + hubbard.ll * yj;
                            const vecXd r = (vecXd(2) << (xi-xj), (yi-yj)).finished();
                            // TODO: check complex or real
                            tmpFourier += cos(-r.dot(q)) * Matsubara(j, i) / hubbard.ls;
                        }
                    }
                }
            }

            obs_mean_gt0_k[l] += tmpFourier;
            obs_err_gt0_k[l] += tmpFourier * tmpFourier;
        }
        obs_mean_gt0_k[l] /= nbin;
        obs_err_gt0_k[l] /= nbin;
        obs_err_gt0_k[l] = pow(obs_err_gt0_k[l]-pow(obs_mean_gt0_k[l], 2), 0.5) / pow(nbin - 1, 0.5);
    }
}
