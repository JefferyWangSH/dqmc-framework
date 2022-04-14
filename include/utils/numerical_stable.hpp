#ifndef UTIL_NUMERICAL_STABLE_HPP
#define UTIL_NUMERICAL_STABLE_HPP
#pragma once

/**
  *  This head file defines the interface Utils::NumericalStable class,
  *  which contains subroutines to help compute equal-time and 
  *  time-displaced (dynamical) Greens function in a stable manner.
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <Eigen/QR>
#include "svd_stack.h"


namespace Utils {

    // ------------------------------- Utils::NumericalStable class ----------------------------------
    // including static subroutines for numerical stablization
    class NumericalStable {
        public:

        using Matrix = Eigen::MatrixXd;
        using Vector = Eigen::VectorXd;

        /*
         *  Subroutine to return the difference between two same-size matrices.
         *  Input: umat, vmat
         *  Output: the maximum difference -> error
         */
        static void matrix_compare_error(const Matrix& umat, const Matrix& vmat, double& error) {
            assert( umat.rows() == vmat.rows() );
            assert( umat.cols() == vmat.cols() );
            assert( umat.rows() == umat.cols() );

            const int ndim = (int)umat.rows();
            double tmp_error = 0.0;
            for (int i = 0; i < ndim; ++i) {
                for (int j = 0; j < ndim; ++j) {
                    tmp_error = std::max(tmp_error, std::abs(umat(i, j) - vmat(i, j)));
                }
            }
            error = tmp_error;
        }

        /*
         *  Subroutine to perform the decomposition of a vector, dvec = dmax * dmin,
         *  to ensure all elements that greater than one are in dmax,
         *  and all elements that less than one are in dmin.
         *  Input: dvec
         *  Output: dmax, dmin
         */
        static void div_dvec_max_min(const Vector& dvec, Vector& dmax, Vector& dmin) {
            assert( dvec.size() == dmax.size() );
            assert( dvec.size() == dmin.size() );

            const int ndim = (int)dvec.size();
            for (int i = 0; i < ndim; ++i) {
                assert( dvec(i) >= 0 );
                if (dvec(i) >= 1.0) {
                    dmin(i) = 1.0; dmax(i) = dvec(i);
                }
                if (dvec(i) < 1.0) {
                    dmax(i) = 1.0; dmin(i) = dvec(i);
                }
            }
        }

        /*
         * Subroutine to perform full matrix * diagonal matrix ^ -1 * full matrix
         * Input: vmat, dvec, umat
         * Output: zmat
         */
        static void mult_v_invd_u(const Matrix& vmat, const Vector& dvec, const Matrix& umat, Matrix& zmat) {
            assert( vmat.cols() == umat.cols() );
            assert( vmat.cols() == zmat.cols() );
            assert( vmat.rows() == umat.rows() );
            assert( vmat.rows() == zmat.rows() );
            assert( vmat.rows() == vmat.cols() );
            assert( vmat.cols() == dvec.size() );

            const int ndim = (int)vmat.rows();

            for (int i = 0; i < ndim; ++i) {
                for (int j = 0; j < ndim; ++j) {
                    double ztmp = 0.0;
                    for (int k = 0; k < ndim; ++k) {
                        ztmp += vmat(j, k) * umat(k, i) / dvec(k);
                    }
                    zmat(j, i) = ztmp;
                }
            }
        }

        /*
         * Subroutine to perform full matrix * diagonal matrix * full matrix
         * Input: vmat, dvec, umat
         * Output: zmat
         */
        static void mult_v_d_u(const Matrix& vmat, const Vector& dvec, const Matrix& umat, Matrix& zmat) {
            assert( vmat.cols() == umat.cols() );
            assert( vmat.cols() == zmat.cols() );
            assert( vmat.rows() == umat.rows() );
            assert( vmat.rows() == zmat.rows() );
            assert( vmat.rows() == vmat.cols() );
            assert( vmat.cols() == dvec.size() );

            const int ndim = (int)vmat.rows();

            for (int i = 0; i < ndim; ++i) {
                for (int j = 0; j < ndim; ++j) {
                    double ztmp = 0.0;
                    for (int k = 0; k < ndim; ++k) {
                        ztmp += vmat(j, k) * umat(k, i) * dvec(k);
                    }
                    zmat(j, i) = ztmp;
                }
            }
        }

        /*
         * return (1 + USV^T)^-1, with method of QR decomposition
         * if dynamical observables measured, return (1 + USV^T)^-1 * USV^T as well.
         */
        static void compute_greens_00_bb(const Matrix& U, const Vector& S, const Matrix& V, Matrix& gtt) {
            // split S = Sbi^-1 * Ss
            Vector Sbi(S.size());
            Vector Ss(S.size());
            for (int i = 0; i < S.size(); ++i) {
                assert( S(i) >= 0 );
                if(S(i) > 1) {
                    Sbi(i) = 1.0/S(i); Ss(i) = 1.0;
                }
                else {
                    Sbi(i) = 1.0; Ss(i) = S(i);
                }
            }

            // compute (1 + USV^T)^-1 in a stable manner
            // note that H is kinda good conditioned, containing small scale information only.
            Matrix H = Sbi.asDiagonal() * U.transpose() + Ss.asDiagonal() * V.transpose();

            /* gtt */
            gtt = H.fullPivHouseholderQr().solve(Sbi.asDiagonal() * U.transpose());
        }

        /*
         * return (1 + USV^T)^-1 * USV^T, with method of QR decomposition
         * to obtain time-displaced Greens function G(\beta, 0)
         */
        static void compute_greens_b0(const Matrix& U, const Vector& S, const Matrix& V, Matrix& gt0) {
            // split S = Sbi^-1 * Ss
            Vector Sbi(S.size());
            Vector Ss(S.size());
            for (int i = 0; i < S.size(); ++i) {
                assert( S(i) >= 0 );
                if(S(i) > 1) {
                    Sbi(i) = 1.0/S(i); Ss(i) = 1.0;
                }
                else {
                    Sbi(i) = 1.0; Ss(i) = S(i);
                }
            }

            // compute (1 + USV^T)^-1 * USV^T in a stable manner
            // note that H is kinda good conditioned, containing small scale information only.
            Matrix H = Sbi.asDiagonal() * U.transpose() + Ss.asDiagonal() * V.transpose();

            /* gt0 */
            gt0 = H.fullPivHouseholderQr().solve(Ss.asDiagonal() * V.transpose());
        }

        /*
         *  return (1 + left * right^T)^-1 in a stable manner, with method of MGS factorization
         *  note: (1 + left * right^T)^-1 = (1 + (USV^T)_left * (VSU^T)_right)^-1
         */
        static void compute_greens_eqtime(SvdStack& left, SvdStack& right, Matrix &gtt) {
            assert(left.MatDim() == right.MatDim());
            const int ndim = left.MatDim();

            /* at l = 0 */
            if ( left.empty() ) {
                compute_greens_00_bb(right.MatrixV(), right.SingularValues(), right.MatrixU(), gtt);
                return;
            }

            /* at l = lt */
            if ( right.empty() ) {
                compute_greens_00_bb(left.MatrixU(), left.SingularValues(), left.MatrixV(), gtt);
                return;
            }

            // local params
            const Matrix ul = left.MatrixU();
            const Vector dl = left.SingularValues();
            const Matrix vl = left.MatrixV();
            const Matrix ur = right.MatrixU();
            const Vector dr = right.SingularValues();
            const Matrix vr = right.MatrixV();

            Vector dlmax(dl.size()), dlmin(dl.size());
            Vector drmax(dr.size()), drmin(dr.size());

            Matrix Atmp(ndim, ndim), Btmp(ndim, ndim);
            Matrix tmp(ndim, ndim);

            /** Modified Gram-Schmidt (MGS) factorization */
            // breakup dr = drmax * drmin , dl = dlmax * dlmin
            div_dvec_max_min(dl, dlmax, dlmin);
            div_dvec_max_min(dr, drmax, drmin);

            // Atmp = ul^T * ur, Btmp = vl^T * vr
            Atmp = ul.transpose() * ur;
            Btmp = vl.transpose() * vr;

            // Atmp = dlmax^-1 * (ul^T * ur) * drmax^-1
            // Btmp = dlmin * (vl^T * vr) * drmin
            for (int j = 0; j < ndim; ++j) {
                for (int i = 0; i < ndim; ++i) {
                    Atmp(i, j) = Atmp(i, j) / (dlmax(i) * drmax(j));
                    Btmp(i, j) = Btmp(i, j) * dlmin(i) * drmin(j);
                }
            }

            tmp = Atmp + Btmp;
            mult_v_invd_u(ur, drmax, tmp.inverse(), Atmp);

            /* gtt */
            mult_v_invd_u(Atmp, dlmax, ul.transpose(), gtt);
        }

        /*
         *  return time-displaced Greens function in a stable manner,
         *  with method of MGS factorization
         */
        static void compute_greens_dynamic(SvdStack& left, SvdStack& right, Matrix &gt0, Matrix &g0t) {
            assert( left.MatDim() == right.MatDim() );
            const int ndim = left.MatDim();

            /* at l = 0 */
            if( left.empty() ) {
                // gt0 = gtt at t = 0
                compute_greens_00_bb(right.MatrixV(), right.SingularValues(), right.MatrixU(), gt0);

                // g0t = - ( 1 - gtt ）at t = 0
                // for convenience: in fact no definition for g0t at t = 0.
                g0t = - (Matrix::Identity(ndim, ndim) - gt0);
                return;
            }

            /* at l = lt */
            if( right.empty() ) {
                // gt0 = ( 1 + B(\beta, 0) )^-1 * B(\beta, 0)
                compute_greens_b0(left.MatrixU(), left.SingularValues(), left.MatrixV(), gt0);

                // g0t = -gtt at t = beta
                compute_greens_00_bb(left.MatrixU(), left.SingularValues(), left.MatrixV(), g0t);
                g0t = - g0t;
                return;
            }

            // local params
            const Matrix ul = left.MatrixU();
            const Vector dl = left.SingularValues();
            const Matrix vl = left.MatrixV();
            const Matrix ur = right.MatrixU();
            const Vector dr = right.SingularValues();
            const Matrix vr = right.MatrixV();

            Vector dlmax(dl.size()), dlmin(dl.size());
            Vector drmax(dr.size()), drmin(dr.size());

            Matrix Atmp(ndim, ndim), Btmp(ndim, ndim);
            Matrix Xtmp(ndim, ndim), Ytmp(ndim, ndim);
            Matrix tmp(ndim, ndim);

            /** Modified Gram-Schmidt (MGS) factorization */
            // breakup dr = drmax * drmin , dl = dlmax * dlmin
            div_dvec_max_min(dl, dlmax, dlmin);
            div_dvec_max_min(dr, drmax, drmin);

            /** gt0 */
            // Atmp = ul^T * ur, Btmp = vl^T * vr
            Atmp = ul.transpose() * ur;
            Btmp = vl.transpose() * vr;

            // Atmp = dlmax^-1 * (ul^T * ur) * drmax^-1
            // Btmp = dlmin * (vl^T * vr) * drmin
            for (int j = 0; j < ndim; ++j) {
                for (int i = 0; i < ndim; ++i) {
                    Atmp(i, j) = Atmp(i, j) / (dlmax(i) * drmax(j));
                    Btmp(i, j) = Btmp(i, j) * dlmin(i) * drmin(j);
                }
            }
            tmp = Atmp + Btmp;
            mult_v_invd_u(ur, drmax, tmp.inverse(), Atmp);
            mult_v_d_u(Atmp, dlmin, vl.transpose(), gt0);

            /** g0t */
            // Xtmp = vr^T * vl, Ytmp = ur^T * ul
            Xtmp = vr.transpose() * vl;
            Ytmp = ur.transpose() * ul;

            // Xtmp = drmax^-1 * (vr^T * vl) * dlmax^-1
            // Ytmp = drmin * (ur^T * ul) * dlmin
            for (int j = 0; j < ndim; ++j) {
                for (int i = 0; i < ndim; ++i) {
                    Xtmp(i, j) = Xtmp(i, j) / (drmax(i) * dlmax(j));
                    Ytmp(i, j) = Ytmp(i, j) * drmin(i) * dlmin(j);
                }
            }
            tmp = Xtmp + Ytmp;
            mult_v_invd_u(-vl, dlmax, tmp.inverse(), Xtmp);
            mult_v_d_u(Xtmp, drmin, ur.transpose(), g0t);
        }


    };

} // namespace Utils

#endif // UTIL_NUMERICAL_STABLE_HPP