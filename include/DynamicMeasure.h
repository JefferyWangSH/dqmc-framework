#ifndef DQMC_HUBBARD_DYNAMICMEASURE_H
#define DQMC_HUBBARD_DYNAMICMEASURE_H
#pragma once

/**
  *  This head file includes module for time-displaced (dynamic) measuring.
  *  Class: measure::dynamicMeasure
  *  Measuring:
  *   1. Dynamical correlation function of imaginary time: G(k, \tua) = < c(k, \tau) * c^+(k, 0) >
  *   2. Helicity modulus \rho_s of superconducting: \rho_s = (\Gamma_L - \Gamma_T) / 4
  *   3. ...
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <vector>

// forward declaration
namespace Model { class Hubbard; }


namespace Measure{

    class DynamicMeasure {
    public:
        int nbin{20};

        /* for time-displaced measurements */

        // dynamical correlation function of imaginary time G(k, \tua) = < c(k, \tau) * c^+(k, 0) >, \tau > 0.
        std::vector<std::vector<double>> bin_g_kt;      // bin_g_kt[bin][tau] <type double>
        std::vector<double> mean_g_kt;                  // mean_g_kt[tau] <type double>
        std::vector<double> err_g_kt;                   // err_g_kt[tau] <type double>

        // helicity modulus \rho_s of superconducting
        std::vector<double> bin_rho_s;                  // bin_rho_s[bin] <type double>
        double mean_rho_s = 0.0;                        // mean_rho_s <type double>
        double err_rho_s = 0.0;                         // err_rho_s <type double>

        std::vector<std::vector<Eigen::MatrixXd>> bin_gt0_up;       // data [bin][tau] <type Eigen::MatrixXd>
        std::vector<std::vector<Eigen::MatrixXd>> bin_g0t_up;
        std::vector<std::vector<Eigen::MatrixXd>> bin_gtt_up;
        std::vector<std::vector<Eigen::MatrixXd>> bin_gt0_dn;
        std::vector<std::vector<Eigen::MatrixXd>> bin_g0t_dn;
        std::vector<std::vector<Eigen::MatrixXd>> bin_gtt_dn;

        // sign problem
        double mean_sign = 0.0;
        double err_sign = 0.0;
        std::vector<double> bin_sign;

        // temporary parameters
        int n_time_displaced = 0;
        double tmp_sign = 0.0;
        std::vector<Eigen::MatrixXd> tmp_gt0_tau_up;        // data[tau] <type Eigen::MatrixXd>
        std::vector<Eigen::MatrixXd> tmp_g0t_tau_up;
        std::vector<Eigen::MatrixXd> tmp_gtt_tau_up;
        std::vector<Eigen::MatrixXd> tmp_gt0_tau_dn;
        std::vector<Eigen::MatrixXd> tmp_g0t_tau_dn;
        std::vector<Eigen::MatrixXd> tmp_gtt_tau_dn;

        // lattice momentum q
        Eigen::VectorXd q = Eigen::VectorXd::Zero(2);

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /* (de)construct functions */
        DynamicMeasure() = default;
        explicit DynamicMeasure(const int &nbin);
        ~DynamicMeasure() = default;

        /* resize the size of bins */
        void resize(const int &nbin);

        /* prepare for measuring */
        void initial(const Model::Hubbard &hubbard);

        /* clear temporary parameters */
        void clear();
        void clear_temporary(const Model::Hubbard &hubbard);

        /* time-displaced measurements */
        void measure_time_displaced(const Model::Hubbard &hubbard);

        /* normalize data from scratch */
        void normalizeStats(const Model::Hubbard &hubbard);

        /* bin measurements */
        void write_Stats_to_bins(const int &bin, const Model::Hubbard &hubbard);

        /* analyse dynamical statistics */
        void analyse_timeDisplaced_Stats(const Model::Hubbard &hubbard);

    private:
        /** dynamical correlation function in momentum space < c(k,tau) * c^+(k,0) > */
        void analyse_Dynamical_Corr(const int &bin, const Model::Hubbard &hubbard);

        /** helicity modules \rho_s */
        void analyse_Rho_S(const int &bin, const Model::Hubbard &hubbard);
    };
}

#endif //DQMC_HUBBARD_DYNAMICMEASURE_H
