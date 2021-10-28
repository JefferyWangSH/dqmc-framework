#ifndef MATRIX_UTIL_HPP
#define MATRIX_UTIL_HPP
#pragma once

/**
 *  This source file includes some diagonalizing tools with C++/Eigen interface
 *  for diagonalizing real matrices using mkl and lapack.
 *  including:
 *    1. generalized SVD decomposition for arbitrary M * N matrices,
 *    2. optimized diagonalizing mechanism for N * N real symmetric matrix
 *  and calculation accuracy is guaranteed.
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
     * @param m -> number of rows.
     * @param n -> number of cols.
     * @param a -> arbitrary M*N real matrix to be solved.
     * @param u -> u matrix in Eigen::Matrix, M * M.
     * @param s -> eigenvalues s in Eigen::Vector, descending sorted.
     * @param v -> v matrix in Eigen::Matrix, N * N.
     */
    void mkl_lapack_dgesvd(const int &m, const int &n, const Eigen::MatrixXd &a, Eigen::MatrixXd &u, Eigen::VectorXd &s, Eigen::MatrixXd &v) {
        assert( m == a.rows() );
        assert( n == a.cols() );

        // Matrix size
        int matrix_layout = LAPACK_ROW_MAJOR;
        lapack_int info, lda = m, ldu = m, ldvt = n;

        // Local arrays
        double _s[ldu * ldu], _u[ldu * m], _vt[ldvt * n];
        double _a[lda * n];
        double superb[ldu * lda];
        for (int i = 0; i < lda * n; ++i) {
            _a[i] = a(i/lda, i%lda);
        }

        // Compute SVD
        info = LAPACKE_dgesvd( matrix_layout, 'A', 'A', m, n, _a, lda, _s, _u, ldu, _vt, ldvt, superb );

        // Check for convergence
        if( info > 0 ) {
            std::cerr << "The algorithm computing SVD failed to converge." << std::endl;
            exit( 1 );
        }

        // Convert results to Eigen
        u = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(_u, n, n);
        s = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(_s, 1, n);
        v = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(_vt, m, m);
    }


    /**
     * Calculate eigenvalues and eigenstates given an arbitrary N * N real symmetric matrix, using MKL_LAPACK
     *      A  ->  T^dagger * S * T
     * where T is rotation matrix, which is orthogonal;
     *       S is diagonal matrix with eigenvalues being diagonal elements.
     *
     * @param N -> number of rows/cols.
     * @param a -> arbitrary N * N real symmetric matrix to be solved.
     * @param s -> diagonal eigen matrix.
     * @param T -> rotation matrix, columns being eigenstates.
     */
    void mkl_lapack_dsyev(const int &N, const Eigen::MatrixXd &a, Eigen::VectorXd &s, Eigen::MatrixXd &T) {
        assert( a.rows() == N );
        assert( a.cols() == N );
        // make sure that a is symmetric
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                assert( a(i, j) == a(j, i) );
            }
        }

        // Locals params
        lapack_int n = N, lda = N, info;
        double w[N];
        double _a[lda * n];
        const Eigen::MatrixXd a_upper = a.triangularView<Eigen::Upper>();

        // Convert eigen matrix to array ( upper triangular )
        for (int i = 0; i < lda * n; ++i) {
            _a[i] = a_upper(i/lda, i%lda);
        }

        // Solve eigen problem
        info = LAPACKE_dsyev( LAPACK_ROW_MAJOR, 'V', 'U', n, _a, lda, w );

        // Check for convergence
        if( info > 0 ) {
            printf( "The algorithm failed to compute eigenvalues.\n" );
            exit( 1 );
        }

        // Convert eigenvalues and eigenvectors to eigen style
        s = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(w, 1, n);
        T = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(_a, n, n);
    }

} // namespce MatrixUtil

#endif //MATRIX_UTIL_HPP