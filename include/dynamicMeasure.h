#ifndef DQMC_HUBBARD_DYNAMICMEASURE_H
#define DQMC_HUBBARD_DYNAMICMEASURE_H
#pragma once

/**
 *  This head file includes module for time-displaced (dynamic) measuring.
 *  Class: measure::dynamicMeasure.
 */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

typedef Eigen::MatrixXd matXd;
typedef Eigen::VectorXd vecXd;

// forward declaration
class Hubbard;

namespace measure{
    class dynamicMeasure {
    public:
        int nbin{20};

        // for time-displaced measurements: Matsubara Green's function gt0
        matXd obs_bin_gt0[500][500];        // obs_bin_gt0 [bin][tau] <type matXd>
        double obs_mean_gt0_k[500]{};       // mean value: obs_mean_gt0_k[tau] <double>
        double obs_err_gt0_k[500]{};        // statistical error: obs_err_gt0_k[tau] <double>

        // temporary parameters
        int n_time_displaced = 0;
        matXd vec_gt0_tau[500];

        // lattice momentum q
        vecXd q = vecXd::Zero(2);

        /* construct function */
        dynamicMeasure() = default;

        /* resize the size of bins */
        void resize(const int &nbin);

        /* prepare for measuring */
        void initial(const Hubbard &hubbard);

        /* clear temporary parameters */
        void clear(const Hubbard &hubbard);

        /* time-displaced measurements */
        void measure_time_displaced(const Hubbard &hubbard);

        /* normalize data from scratch */
        void normalizeStats(const Hubbard &hubbard);

        /* bin measurements */
        void write_Stats_to_bins(const int &bin, const Hubbard &hubbard);

        /** time-displaced green function in momentum space < c(k,tau) * c^T(k,0) > */
        void analyse_timeDisplaced_Stats(const Hubbard &hubbard);
    };
}

#endif //DQMC_HUBBARD_DYNAMICMEASURE_H
