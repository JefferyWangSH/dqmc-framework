#ifndef HUBBARD_V1_3_HUBBARD_H
#define HUBBARD_V1_3_HUBBARD_H
#pragma once

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
    int ll{4}, ls{16}, lt{40};
    double beta{4.0}, dtau{0.1};
    double t{1.0}, Uint{4.0}, mu{0.0}, alpha{};

    int nwrap{8};
    int current_tau{0};

    matXd s, expmdtK, exppdtK;      // aux field and kinetic matrix expK
    matXd GreenU, GreenD;
    std::vector<matXd> vecGreenU, vecGreenD;

    // aux SvdStack class for stable multiplication
    SvdStack stackLeftU;
    SvdStack stackLeftD;
    SvdStack stackRightU;
    SvdStack stackRightD;

public:
    /** construction */
    Hubbard() = default;
    Hubbard(int ll, int lt, double beta,
            double t, double Uint, double mu, int nwrap=8);

    /** sweep the space-time lattice from 0 to beta */
    void sweep_0_to_beta(int istab);

    /** sweep the space-time lattice from beta to 0 */
    void sweep_beta_to_0(int istab);

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

#endif //HUBBARD_V1_3_HUBBARD_H
