#ifndef DQMC_HUBBARD_CHECKERBOARD_H
#define DQMC_HUBBARD_CHECKERBOARD_H
#pragma once

/**
 *  This header file includes class and subroutines for multiplying hopping matrix expK to a dense matrix,
 *  with high computational efficiency by method of checkerboard break-up for even lattice sizes,
 *  and by direct multiplication for others.
 */


#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

// forward declaration
class Hubbard;

class CheckerBoard {

public:
    // lattice params
    int ll{}, ls{};
    double t{}, mu{};
    double dtau{};

    // helper matrices for constructing exponent of hopping K
    Eigen::MatrixXd exp_dtK, inv_exp_dtK;
    Eigen::Matrix4d exp_dtK_reduced, inv_exp_dtK_reduced;

    bool is_checkerboard = true;

public:
    /** construction */
    CheckerBoard() = default;

    /** return is_checkerboard */
    bool is_checker_board() const;

    /** initialize from Hubbard class */
    void init_from_model(const Hubbard &hubbard);

    /** high-efficient multiplication of hopping K */
    void mult_expK_from_left(Eigen::MatrixXd &A) const;

    void mult_inv_expK_from_left(Eigen::MatrixXd &A) const;

    void mult_expK_from_right(Eigen::MatrixXd &A) const;

    void mult_inv_expK_from_right(Eigen::MatrixXd &A) const;

    void mult_trans_expK_from_left(Eigen::MatrixXd &A) const;

private:

    /** checkout prerequisite */
    void check_checker_board();

    /** generate exponent of hopping K by direct multiplication */
    void make_expK_direct();

    /** generate exponent of hopping K by checkerboard break-ups */
    void make_expK_checkerboard();

    /** multiplication of hopping K within single plaquette, labeled by site (x,y) */
    void mult_expK_plaquette_from_left(Eigen::MatrixXd &A, int x, int y) const;

    void mult_inv_expK_plaquette_from_left(Eigen::MatrixXd &A, int x, int y) const;

    void mult_expK_plaquette_from_right(Eigen::MatrixXd &A, int x, int y) const;

    void mult_inv_expK_plaquette_from_right(Eigen::MatrixXd &A, int x, int y) const;
};

#endif //DQMC_HUBBARD_CHECKERBOARD_H
