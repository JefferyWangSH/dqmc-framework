#ifndef DQMC_HUBBARD_EQTIMEMEASURE_H
#define DQMC_HUBBARD_EQTIMEMEASURE_H
#pragma once

/**
  *  This head file includes module for equal-time measuring.
  *  Class: Measure::eqtimeMeasure
  *  Measuring:
  *   1. Double occupancy D = < n_up*n_dn >
  *   2. Single particle kinetic energy
  *   3. Distributions of electrons in momentum space
  *   4. Local spin correlation, magnetization C(0,0) = < (n_up - n_dn)^2 >
  *   5. Spin density structure factor (SDW)
  *   6. Charge density structure factor (CDW)
  *   6. Space correlation of s-wave pairing
  *   7. ...
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <vector>
#include "measure_data.h"

// forward declaration
namespace Model { class Hubbard; }
namespace Measure { class MeasureData; }


namespace Measure{

    class EqtimeMeasure {
    public:
        int nbin{20};

        /* for equal-time (static) measurements */
        Measure::MeasureData sign;                              // average sign to keep track of sign problem
        Measure::MeasureData filling_number;                    // filling number <n>
        Measure::MeasureData double_occupancy;                  // double occupancy
        Measure::MeasureData kinetic_energy;                    // kinetic energy
        Measure::MeasureData momentum_distribution;             // particle distribution in momentum space
        Measure::MeasureData local_spin_corr;                   // local spin correlation
        Measure::MeasureData spin_density_structure_factor;     // order parameter of spin density wave (SDW)
        Measure::MeasureData charge_density_structure_factor;   // order parameter of charge density wave (CDW)
        std::vector<Measure::MeasureData> pairing_corr;         // space correlation of s-wave pairing

        // lattice momentum q
        Eigen::VectorXd q = Eigen::VectorXd::Zero(2);

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /* (de)construct functions */
        EqtimeMeasure()  = default;
        explicit EqtimeMeasure(const int &nbin);

        ~EqtimeMeasure() = default;

        /* resize the size of bins */
        void resize(const int &nbin);

        /* prepare for measuring */
        void initial(const Model::Hubbard &hubbard);

        /* clear temporary parameters */
        void clear_temporary();

        /* equal-time measurements */
        void equal_time_measure(const Model::Hubbard &hubbard);

        /* normalize data from scratch */
        void normalize_stats(const Model::Hubbard &hubbard);

        /* bin measurements */
        void write_stats_to_bins(int bin);

        /* analyse all equal-time statistics from simulation */
        void analyse_stats(const Model::Hubbard &hubbard);

    private:
        /* filling number: <n> = \sum_i ( n_i_up + n_i_dn ) */
        void measure_filling_number(const int &tau, const Model::Hubbard &hubbard);

        /** double occupation: D = < n_up * n_dn > */
        void measure_double_occupancy(const int &tau, const Model::Hubbard &hubbard);

        /** single particle kinetic energy */
        void measure_kinetic_energy(const int &tau, const Model::Hubbard &hubbard);

        /** momentum distribution of electrons: fourier transformation of real-space electron distribution */
        void measure_momentum_distribution(const int &tau, const Model::Hubbard &hubbard);

        /** local spin correlation: magnetization C(0,0) = < (n_up - n_dn)^2 > */
        void measure_local_spin_corr(const int &tau, const Model::Hubbard &hubbard);

        /** spin density wave (SDW) structure factor S(q): fourier transformation of real-space correlation of spins */
        void measure_spin_density_structure_factor(const int &tau, const Model::Hubbard &hubbard);

        /* charge density wave (CDW) structure factor C(q) */
        void measure_charge_density_structure_factor(const int &tau, const Model::Hubbard &hubbard);

        /** space correlation of s-wave pairing order parameter */
        void measure_pairing_corr(const int &tau, const Model::Hubbard &hubbard);
    };
}

#endif //DQMC_HUBBARD_EQTIMEMEASURE_H
