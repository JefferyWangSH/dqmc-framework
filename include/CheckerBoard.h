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
    int ll{}, ls{};
    double t{}, mu{};
    double dtau{};

    Eigen::MatrixXd exp_dtK, inv_exp_dtK;
    Eigen::Matrix4d exp_dtK_reduced, inv_exp_dtK_reduced;

    bool is_checkerboard = true;

public:
    CheckerBoard() = default;

    bool is_checker_board() const;

    void init_from_model(const Hubbard &hubbard);

    void mult_expK_from_left(Eigen::MatrixXd &A) const;

    void mult_inv_expK_from_left(Eigen::MatrixXd &A) const;

    void mult_expK_from_right(Eigen::MatrixXd &A) const;

    void mult_inv_expK_from_right(Eigen::MatrixXd &A) const;

    void mult_trans_expK_from_left(Eigen::MatrixXd &A) const;

private:

    void check_checker_board();

    void make_expK_direct();

    void make_expK_checkerboard();

    void mult_expK_plaquette_from_left(Eigen::MatrixXd &A, int x, int y) const;

    void mult_inv_expK_plaquette_from_left(Eigen::MatrixXd &A, int x, int y) const;

    void mult_expK_plaquette_from_right(Eigen::MatrixXd &A, int x, int y) const;

    void mult_inv_expK_plaquette_from_right(Eigen::MatrixXd &A, int x, int y) const;
};

#endif //DQMC_HUBBARD_CHECKERBOARD_H
