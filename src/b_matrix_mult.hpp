#ifndef DQMC_HUBBARD_BMATRIXMULT_HPP
#define DQMC_HUBBARD_BMATRIXMULT_HPP
#pragma once

/**
  *  This source file includes subroutines for multiplication of
  *  B matrices during configuration updates and green function fabrication.
  *  Both direct multiplication and checkerboard break-up
  *  algorithm (for even lattice sizes) are supported.
  */

#include "hubbard.h"

//Eigen::MatrixXd Model::Hubbard::make_B_l(int l, int sigma) {
//    /*
//     *  Compute B matrix at time slice l for given spin-1/2 state.
//     *  definition of B(l, sigma):
//     *  B_l = exp(-\\Delta \\tau V^{\\sigma}(l)) * exp(-\\Delta \\tau K)
//     */
//    assert( l >= 0 && l <= lt );
//    assert( sigma == 1 || sigma == -1 );
//
//    const int tau = (l==0)? lt-1 : l-1;
//    const int eff_sigma = (is_attractive_u)? +1 : sigma;
//
//    Eigen::MatrixXd r;
//    if (this->checkerboard->is_checker_board()) {
//        r = Eigen::MatrixXd::Identity(ls, ls);
//        this->checkerboard->mult_expK_from_left(r);
//    }
//    else { r = this->checkerboard->exp_dtK;}
//
//    for (int i = 0; i < ls; ++i) {
//        r.row(i) *= exp(+ eff_sigma * this->alpha * (*this->s)(i, tau));
//    }
//    return r;
//}

void Model::Hubbard::mult_B_from_left(Eigen::MatrixXd &A, int l, int sigma) {
    /*
     *  Multiply a dense matrix A from the left by B_l
     *  A ->  B_l * A = exp(-\\Delta \\tau V^{\\sigma}(l)) * exp(-\\Delta \\tau K) * A
     *  Matrix A is changed in place.
     */
    assert( A.rows() == this->ls && A.cols() == this->ls );
    assert( l >= 0 && l <= this->lt );
    assert( sigma == 1 || sigma == -1 );

    const int tau = (l==0)? this->lt-1 : l-1;
    const int eff_sigma = (this->is_attractive_u)? +1 : sigma;
    this->checkerboard->mult_expK_from_left(A);
    for (int i = 0; i < this->ls; ++i) {
        A.row(i) *= exp(+ eff_sigma * this->alpha * (*this->s)(i, tau));
    }
}

void Model::Hubbard::mult_B_from_right(Eigen::MatrixXd &A, int l, int sigma) {
    /*
     *  Multiply a dense matrix A from the right by B_l
     *  A ->  A * B_l = A * exp(-\\Delta \\tau V^{\\sigma}(l)) * exp(-\\Delta \\tau K)
     *  Matrix A is changed in place.
     */
    assert( A.rows() == this->ls && A.cols() == this->ls );
    assert( l >= 0 && l <= this->lt );
    assert( sigma == 1 || sigma == -1 );

    const int tau = (l==0)? this->lt-1 : l-1;
    const int eff_sigma = (this->is_attractive_u)? +1 : sigma;
    for (int i = 0; i < this->ls; ++i) {
        A.col(i) *= exp(+ eff_sigma * this->alpha * (*this->s)(i, tau));
    }
    this->checkerboard->mult_expK_from_right(A);
}

void Model::Hubbard::mult_invB_from_left(Eigen::MatrixXd &A, int l, int sigma) {
    /*
     *  Multiply a dense matrix A from the left by B_l^{-1}
     *  A -> B_{l}^{-1} * A = exp(+\\Delta \\tau K) * exp(+\\Delta \\tau V^{\\sigma}(l))  * A
     *  Matrix A is changed in place.
     */
    assert( A.rows() == this->ls && A.cols() == this->ls );
    assert( l >= 0 && l <= this->lt );
    assert( sigma == 1 || sigma == -1 );

    const int tau = (l==0)? this->lt-1 : l-1;
    const int eff_sigma = (this->is_attractive_u)? +1 : sigma;
    for (int i = 0; i < this->ls; ++i) {
        A.row(i) *= exp(- eff_sigma * this->alpha * (*this->s)(i, tau));
    }
    this->checkerboard->mult_inv_expK_from_left(A);
}

void Model::Hubbard::mult_invB_from_right(Eigen::MatrixXd &A, int l, int sigma) {
    /*
     *  Multiply a dense matrix A from the right by B_l^{-1}
     *  A -> A * B_{l}^{-1} = A * exp(+\\Delta \\tau K) * exp(+\\Delta \\tau V^{\\sigma}(l))
     *  Matrix A is changed in place.
     */
    assert( A.rows() == this->ls && A.cols() == this->ls );
    assert( l >= 0 && l <= this->lt );
    assert( sigma == 1 || sigma == -1 );

    const int tau = (l==0)? this->lt-1 : l-1;
    const int eff_sigma = (this->is_attractive_u)? +1 : sigma;
    this->checkerboard->mult_inv_expK_from_right(A);
    for (int i = 0; i < this->ls; ++i) {
        A.col(i) *= exp(- eff_sigma * this->alpha * (*this->s)(i, tau));
    }
}

void Model::Hubbard::mult_transB_from_left(Eigen::MatrixXd &A, int l, int sigma) {
    /*
     *  Multiply a dense matrix A from the left by B_l^T
     *  A ->  B_l^T * A = exp(-\\Delta \\tau K)^T * exp(-\\Delta \\tau V^{\\sigma}(l)) * A
     *  Matrix A is changed in place.
     */
    assert( A.rows() == this->ls && A.cols() == this->ls );
    assert( l >= 0 && l <= this->lt );
    assert( sigma == 1 || sigma == -1 );

    const int tau = (l==0)? this->lt-1 : l-1;
    const int eff_sigma = (this->is_attractive_u)? +1 : sigma;
    for (int i = 0; i < this->ls; ++i) {
        A.row(i) *= exp(+ eff_sigma * this->alpha * (*this->s)(i, tau));
    }
    this->checkerboard->mult_trans_expK_from_left(A);
}

#endif //DQMC_HUBBARD_BMATRIXMULT_HPP