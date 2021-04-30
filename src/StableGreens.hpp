#ifndef DQMC_HUBBARD_STABLEGREENS_HPP
#define DQMC_HUBBARD_STABLEGREENS_HPP
#pragma once

/**
 *  This head file includes subroutines
 *  to help compute equal-time Greens function in a stable manner.
 */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include "SvdStack.hpp"


void matrix_compare_error(const Eigen::MatrixXd &umat, const Eigen::MatrixXd &vmat, double &error) {
    /*
     *  Subroutine to return the difference between two same-size matrices.
     *  Input: umat, vmat
     *  Output: the maximum difference -> error
     */
    assert(umat.rows() == vmat.rows());
    assert(umat.cols() == vmat.cols());
    assert(umat.rows() == umat.cols());

    const int ndim = (int)umat.rows();
    double tmp_error = 0.0;
    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) {
            tmp_error = std::max(tmp_error, std::abs(umat(i, j)-vmat(i, j)));
        }
    }
    error = tmp_error;
}

void div_dvec_max_min(const Eigen::VectorXd &dvec, Eigen::VectorXd &dmax, Eigen::VectorXd &dmin) {
    /*
     *  Subroutine to perform the decomposition of a vector, dvec = dmax * dmin,
     *  to ensure all elements that greater than one are in dmax,
     *  and all elements that less than one are in dmin.
     *  Input: dvec
     *  Output: dmax, dmin
     */

    assert(dvec.size() == dmax.size());
    assert(dvec.size() == dmin.size());

    const int ndim = (int)dvec.size();
    for (int i = 0; i < ndim; ++i) {
        assert(dvec(i) >= 0);
        if (dvec(i) >= 1.0) {
            dmin(i) = 1.0; dmax(i) = dvec(i);
        }
        if (dvec(i) < 1.0) {
            dmax(i) = 1.0; dmin(i) = dvec(i);
        }
    }
}

void mult_v_invd_u(const Eigen::MatrixXd &vmat, const Eigen::VectorXd &dvec, const Eigen::MatrixXd &umat, Eigen::MatrixXd &zmat) {
    /*
     * Subroutine to perform full matrix * diagonal matrix ^ -1 * full matrix
     * Input: vmat, dvec, umat
     * Output: zmat
     */
    assert(vmat.cols() == umat.cols());
    assert(vmat.cols() == zmat.cols());
    assert(vmat.rows() == umat.rows());
    assert(vmat.rows() == zmat.rows());
    assert(vmat.rows() == vmat.cols());
    assert(vmat.cols() == dvec.size());

    const int ndim = (int)vmat.rows();

    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) {
            double ztmp = 0.0;
            for (int k = 0; k < ndim; ++k) {
                ztmp += vmat(j,k) * umat(k,i) / dvec(k);
            }
            zmat(j,i) = ztmp;
        }
    }
}

void mult_v_d_u(const Eigen::MatrixXd &vmat, const Eigen::VectorXd &dvec, const Eigen::MatrixXd &umat, Eigen::MatrixXd &zmat) {
    /*
     * Subroutine to perform full matrix * diagonal matrix * full matrix
     * Input: vmat, dvec, umat
     * Output: zmat
     */
    assert(vmat.cols() == umat.cols());
    assert(vmat.cols() == zmat.cols());
    assert(vmat.rows() == umat.rows());
    assert(vmat.rows() == zmat.rows());
    assert(vmat.rows() == vmat.cols());
    assert(vmat.cols() == dvec.size());

    const int ndim = (int)vmat.rows();

    for (int i = 0; i < ndim; ++i) {
        for (int j = 0; j < ndim; ++j) {
            double ztmp = 0.0;
            for (int k = 0; k < ndim; ++k) {
                ztmp += vmat(j,k) * umat(k,i) * dvec(k);
            }
            zmat(j,i) = ztmp;
        }
    }
}

void compute_Green_00_bb(const Eigen::MatrixXd &U, const Eigen::VectorXd &S, const Eigen::MatrixXd &V, Eigen::MatrixXd &gtt) {
    /*
     * returns (1 + USV^T)^-1, with method of QR decomposition
     * if bool_measure_dynamic, return (1 + USV^T)^-1 * USV^T too.
     */

    // split S = Sbi^-1 * Ss
    Eigen::VectorXd Sbi(S.size());
    Eigen::VectorXd Ss(S.size());
    for (int i = 0; i < S.size(); ++i) {
        assert(S(i) >= 0);
        if(S(i) > 1) {
            Sbi(i) = 1.0/S(i); Ss(i) = 1.0;
        }
        else {
            Sbi(i) = 1.0; Ss(i) = S(i);
        }
    }

    // compute (1 + USV^T)^-1 in a stable manner
    // note that H is kinda good conditioned, containing small scale information only.
    Eigen::MatrixXd H = Sbi.asDiagonal() * U.transpose() + Ss.asDiagonal() * V.transpose();

    /* gtt */
    gtt = H.fullPivHouseholderQr().solve(Sbi.asDiagonal() * U.transpose());
}

void compute_Green_b0(const Eigen::MatrixXd &U, const Eigen::VectorXd &S, const Eigen::MatrixXd &V, Eigen::MatrixXd &gt0) {
    /*
     * returns (1 + USV^T)^-1 * USV^T, with method of QR decomposition
     * to obtain time-displaced Greens function G(\beta, 0)
     */

    // split S = Sbi^-1 * Ss
    Eigen::VectorXd Sbi(S.size());
    Eigen::VectorXd Ss(S.size());
    for (int i = 0; i < S.size(); ++i) {
        assert(S(i) >= 0);
        if(S(i) > 1) {
            Sbi(i) = 1.0/S(i); Ss(i) = 1.0;
        }
        else {
            Sbi(i) = 1.0; Ss(i) = S(i);
        }
    }

    // compute (1 + USV^T)^-1 * USV^T in a stable manner
    // note that H is kinda good conditioned, containing small scale information only.
    Eigen::MatrixXd H = Sbi.asDiagonal() * U.transpose() + Ss.asDiagonal() * V.transpose();

    /* gt0 */
    gt0 = H.fullPivHouseholderQr().solve(Ss.asDiagonal() * V.transpose());
}

void compute_Green_eqtime(const SvdStack *left, const SvdStack *right, Eigen::MatrixXd &gtt) {
    /*
     *  returns (1 + left * right^T)^-1 in a stable manner, with method of MGS factorization
     *  note: (1 + left * right^T)^-1 = (1 + (USV^T)_left * (VSU^T)_right)^-1
     */
    assert(left->n == right->n);
    const int ndim = left->n;

    /* at l = 0 */
    if (left->empty()) {
        compute_Green_00_bb(right->matrixV(), right->singularValues(), right->matrixU(), gtt);
        return;
    }

    /* at l = lt */
    if (right->empty()) {
        compute_Green_00_bb(left->matrixU(), left->singularValues(), left->matrixV(), gtt);
        return;
    }

    // local params
    const Eigen::MatrixXd& ul = left->matrixU();
    const Eigen::VectorXd& dl = left->singularValues();
    const Eigen::MatrixXd& vl = left->matrixV();
    const Eigen::MatrixXd& ur = right->matrixU();
    const Eigen::VectorXd& dr = right->singularValues();
    const Eigen::MatrixXd& vr = right->matrixV();

    Eigen::VectorXd dlmax(dl.size()), dlmin(dl.size());
    Eigen::VectorXd drmax(dr.size()), drmin(dr.size());

    Eigen::MatrixXd Atmp(ndim, ndim), Btmp(ndim, ndim);
    Eigen::MatrixXd tmp(ndim, ndim);

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
            Atmp(i,j) = Atmp(i,j) / (dlmax(i) * drmax(j));
            Btmp(i,j) = Btmp(i,j) * dlmin(i) * drmin(j);
        }
    }

    tmp = Atmp + Btmp;
    mult_v_invd_u(ur, drmax, tmp.inverse(), Atmp);

    /* gtt */
    mult_v_invd_u(Atmp, dlmax, ul.transpose(), gtt);
}

void compute_Green_displaced(const SvdStack *left, const SvdStack *right, Eigen::MatrixXd &gt0, Eigen::MatrixXd &g0t) {
    /*
     *  returns time-displaced Greens function in a stable manner,
     *  with method of MGS factorization
     */
    assert(left->n == right->n);
    const int ndim = left->n;

    /* at l = 0 */
    if(left->empty()) {
        // gt0 = gtt at t = 0
        compute_Green_00_bb(right->matrixV(), right->singularValues(), right->matrixU(), gt0);

        // g0t = - ( 1 - gtt ）at t = 0
        // for convenience: in fact no definition for g0t at t = 0.
        g0t = - (Eigen::MatrixXd::Identity(ndim, ndim) - gt0);
        return;
    }

    /* at l = lt */
    if(right->empty()) {
        // gt0 = ( 1 + B(\beta, 0) )^-1 * B(\beta, 0)
        compute_Green_b0(left->matrixU(), left->singularValues(), left->matrixV(), gt0);

        // g0t = -gtt at t = beta
        compute_Green_00_bb(left->matrixU(), left->singularValues(), left->matrixV(), g0t);
        g0t = -g0t;
        return;
    }

    // local params
    const Eigen::MatrixXd& ul = left->matrixU();
    const Eigen::VectorXd& dl = left->singularValues();
    const Eigen::MatrixXd& vl = left->matrixV();
    const Eigen::MatrixXd& ur = right->matrixU();
    const Eigen::VectorXd& dr = right->singularValues();
    const Eigen::MatrixXd& vr = right->matrixV();

    Eigen::VectorXd dlmax(dl.size()), dlmin(dl.size());
    Eigen::VectorXd drmax(dr.size()), drmin(dr.size());

    Eigen::MatrixXd Atmp(ndim, ndim), Btmp(ndim, ndim);
    Eigen::MatrixXd Xtmp(ndim, ndim), Ytmp(ndim, ndim);
    Eigen::MatrixXd tmp(ndim, ndim);

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
            Atmp(i,j) = Atmp(i,j) / (dlmax(i) * drmax(j));
            Btmp(i,j) = Btmp(i,j) * dlmin(i) * drmin(j);
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
            Xtmp(i,j) = Xtmp(i,j) / (drmax(i) * dlmax(j));
            Ytmp(i,j) = Ytmp(i,j) * drmin(i) * dlmin(j);
        }
    }
    tmp = Xtmp + Ytmp;
    mult_v_invd_u(-vl, dlmax, tmp.inverse(), Xtmp);
    mult_v_d_u(Xtmp, drmin, ur.transpose(), g0t);
}


#endif //DQMC_HUBBARD_STABLEGREENS_HPP