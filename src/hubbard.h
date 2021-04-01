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

#include "svdstack.h"

typedef Eigen::MatrixXd matXd;
typedef Eigen::VectorXd vecXd;

// random engine
static std::default_random_engine gen(time(nullptr));


class Hubbard {

private:
    // model params
    int ll{4}, ls{16}, lt{80};
    double beta{4.0}, dtau{0.1};
    double t{1.0}, Uint{4.0}, mu{0.0}, alpha{};

    bool u_is_attractive = true;

    int nwrap{10};
    int current_tau{0};

    // aux field and kinetic matrix expK
    matXd s, expmdtK, exppdtK;      

    // green function for both spin-1/2 states
    matXd GreenU, GreenD;
    std::vector<matXd> vecGreenU, vecGreenD;

    // time-displaced green function for dynamic measurements
    matXd Green_t0_up, Green_t0_dn; 
    std::vector<matXd> vecGreen_t0_up, vecGreen_t0_dn;

    // aux SvdStack class for stable multiplication
    SvdStack stackLeftU;
    SvdStack stackLeftD;
    SvdStack stackRightU;
    SvdStack stackRightD;

public:
    /** construction */
    Hubbard() = default;
    Hubbard(int ll, int lt, double beta, double t, double Uint, double mu, int nwrap);

    /** sweep the space-time lattice from 0 to beta */
    void sweep_0_to_beta(int istab);

    /** sweep the space-time lattice from beta to 0 */
    void sweep_beta_to_0(int istab);

    /** sweep from beta to 0 to calculate time-displaced green functions */
    void sweep_0_to_beta_displaced(int istab);

    // unnecessary
    friend class detQMC;


private:
    /** randomly initialize aux field */
    void initRandom();

    /** compute exp of Kinetic matrix K */
    void make_expdtK();

    /** initialize udv stacks for sweep use */
    void initStacks(int istab);

    /** compute B matrix with slice l ans spin sigma given*/
    matXd make_Bl(int l, int sigma);

    /** multiply B matrix in place */
    void multB_fromL(matXd& A, int l, int sigma);

    void multB_fromR(matXd& A, int l, int sigma);

    void multinvB_fromL(matXd& A, int l, int sigma);

    void multinvB_fromR(matXd& A, int l, int sigma);

    /** update the aux field at time slice l with Metropolis algorithm */
    void Metropolis_update(int l);

    /** propagate the green's function from l to l+1 */
    void wrap_north(int l);

    /** propagate the green's function from l to l-1 */
    void wrap_south(int l);
};

#endif //DQMC_HUBBARD_HUBBARD_H
