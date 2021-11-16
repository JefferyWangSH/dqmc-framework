#ifndef MATRIX_UTIL_HPP
#define MATRIX_UTIL_HPP
#pragma once

/**
 *  This source file includes some diagonalizing tools with C++/Eigen interface
 *  for diagonalizing real matrices using mkl and lapack.
 *  including:
 *    1. generalized SVD decomposition for arbitrary M * N matrices
 *    2. optimized diagonalizing mechanism for N * N real symmetric matrix
 *  The calculation accuracy and efficiency are guaranteed.
 */


#include <iostream>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "mkl_lapacke.h"

namespace MatrixUtil {

    /**
      * SVD decomposition of arbitrary M * N real matrix, using MKL_LAPACK:
      *      A  ->  U * S * V^T
      * Remind that V is returned in this subroutine, not V transpose.
      *
      * @param row -> number of rows.
      * @param col -> number of cols.
      * @param mat -> arbitrary `row` * `col` real matrix to be solved.
      * @param u -> u matrix in Eigen::Matrix, `row` * `row`.
      * @param s -> eigenvalues s in Eigen::Vector, descending sorted.
      * @param v -> v matrix in Eigen::Matrix, `col` * `col`.
      */
    void mkl_lapack_dgesvd(const int &row, const int &col, const Eigen::MatrixXd &mat, Eigen::MatrixXd &u, Eigen::VectorXd &s, Eigen::MatrixXd &v) {
        assert( row == mat.rows() );
        assert( col == mat.cols() );

        // matrix size
        int matrix_layout = LAPACK_ROW_MAJOR;
        lapack_int info, lda = row, ldu = row, ldvt = col;

        // local arrays
        double _s[ldu * ldu], _u[ldu * row], _vt[ldvt * col];
        double a[lda * col];
        double superb[ldu * lda];
        for (int i = 0; i < lda * col; ++i) {
            a[i] = mat(i/lda, i%lda);
        }

        // compute SVD
        info = LAPACKE_dgesvd( matrix_layout, 'A', 'A', row, col, a, lda, _s, _u, ldu, _vt, ldvt, superb );

        // check for convergence
        if( info > 0 ) {
            std::cerr << "The algorithm computing SVD failed to converge." << std::endl;
            exit( 1 );
        }

        // convert results to Eigen
        u = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(_u, col, col);
        s = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(_s, 1, col);
        v = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(_vt, row, row);
    }


    /**
      * Calculate eigenvalues and eigenstates given an arbitrary N * N real symmetric matrix, using MKL_LAPACK
      *      A  ->  T^dagger * S * T
      * where T is rotation matrix, which is orthogonal;
      *       S is diagonal matrix with eigenvalues being diagonal elements.
      *
      * @param size -> number of rows/cols.
      * @param mat -> arbitrary `size` * `size` real symmetric matrix to be solved.
      * @param s -> diagonal eigen matrix.
      * @param t -> rotation matrix, columns being eigenstates.
      */
    void mkl_lapack_dsyev(const int &size, const Eigen::MatrixXd &mat, Eigen::VectorXd &s, Eigen::MatrixXd &t) {
        assert( mat.rows() == size );
        assert( mat.cols() == size );
        // make sure the input matrix is symmetric
        assert( mat.isApprox(mat.transpose(), 1e-12) );

        // locals params
        lapack_int n = size, lda = size, info;
        double w[n];
        double a[lda * n];

        // convert eigen matrix to c-style array
        for (int i = 0; i < lda * n; ++i) {
            a[i] = mat(i/lda, i%lda);
        }

        // solve eigen problem
        info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, a, lda, w );

        // check for convergence
        if( info > 0 ) {
            std::cerr << "The algorithm failed to compute eigenvalues." << std::endl;
            exit( 1 );
        }

        // convert eigenvalues and eigenvectors to eigen style
        s = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(w, 1, n);
        t = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(a, n, n);
    }

} // namespce MatrixUtil

#endif //MATRIX_UTIL_HPP