#ifndef DQMC_HUBBARD_HUBBARD_H
#define DQMC_HUBBARD_HUBBARD_H
#pragma once

/**
 *  This head file includes hubbard class
 *  which is defined for the DQMC simulation
 *  of half-filled Hubbard model with repulsive interaction potential.
 *  including: model parameters and Monte Carlo update algorithm.
 */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2

#include <Eigen/LU>
#include <Eigen/SVD>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <unsupported/Eigen/MatrixFunctions>
#include <cassert>

#include "SvdStack.h"
#include "eqtimeMeasure.h"
#include "dynamicMeasure.h"


// random engine
static std::default_random_engine gen(time(nullptr));

class Hubbard {

public:
    // model params
    int ll{4}, ls{16}, lt{80};
    double beta{4.0}, dtau{0.1};
    double t{1.0}, Uint{4.0}, mu{0.0}, alpha{};

    bool u_is_attractive = true;

    int nwrap{10};
    int current_tau{0};

    double config_sign{0.0};

    double max_wrap_error_equal{0.0};
    double max_wrap_error_displaced{0.0};

    // aux field and kinetic matrix expK
    Eigen::MatrixXd s, expmdtK, exppdtK;

    // equal-time greens function for both spin-1/2 states
    // critical quantities in DQMC simulation
    Eigen::MatrixXd green_tt_up, green_tt_dn;
    std::vector<Eigen::MatrixXd> vec_green_tt_up, vec_green_tt_dn;

    // time-displaced greens function for dynamic measurements
    // Matsubara greens function: G(\tau, 0) and G(0, \tau).
    // Gij(\tau, 0) = < ci(\tau) * cj^+ (0) >
    Eigen::MatrixXd green_t0_up, green_t0_dn;
    std::vector<Eigen::MatrixXd> vec_green_t0_up, vec_green_t0_dn;

    // Gij(0, \tau) = - < cj^+(\tau) * ci(0) >
    Eigen::MatrixXd green_0t_up, green_0t_dn;
    std::vector<Eigen::MatrixXd> vec_green_0t_up, vec_green_0t_dn;

    // aux SvdStack class for stable multiplication
    SvdStack *stackLeftU{};
    SvdStack *stackLeftD{};
    SvdStack *stackRightU{};
    SvdStack *stackRightD{};

public:
    /** construction */
    Hubbard() = default;
    Hubbard(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap);

    /** sweep the space-time lattice from 0 to beta */
    void sweep_0_to_beta(int is_stable);

    /** sweep the space-time lattice from beta to 0 */
    void sweep_beta_to_0(int is_stable);

    /** sweep from beta to 0 to calculate time-displaced green functions */
    void sweep_0_to_beta_displaced(int is_stable);

    // friend class
    friend class detQMC;
    friend class measure::eqtimeMeasure;
    friend class measure::dynamicMeasure;


private:
    /** randomly initialize aux field */
    void initRandom();

    /** compute exp of Kinetic matrix K */
    void make_expdtK();

    /** initialize udv stacks for sweep use */
    void initStacks(int is_stable);

    /** compute B matrix with slice l ans spin sigma given*/
    Eigen::MatrixXd make_Bl(int l, int sigma);

    /** multiply B matrix in place */
    void multB_fromL(Eigen::MatrixXd& A, int l, int sigma);

    void multB_fromR(Eigen::MatrixXd& A, int l, int sigma);

    void multinvB_fromL(Eigen::MatrixXd& A, int l, int sigma);

    void multinvB_fromR(Eigen::MatrixXd& A, int l, int sigma);

    /** update the aux field at time slice l with Metropolis algorithm */
    void Metropolis_update(int l);

    /** propagate the green's function from l to l+1 */
    void wrap_north(int l);

    /** propagate the green's function from l to l-1 */
    void wrap_south(int l);
};

#endif //DQMC_HUBBARD_HUBBARD_H
