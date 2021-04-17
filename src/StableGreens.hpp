#ifndef DQMC_HUBBARD_STABLEGREENS_HPP
#define DQMC_HUBBARD_STABLEGREENS_HPP
#pragma once

/**
 *  This head file includes subroutines
 *  to help compute equal-time Greens function in a stable manner.
 */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include "svdstack.hpp"

void div_dvec_max_min(const vecXd& dvec, vecXd& dmax, vecXd& dmin) {
    /*
     *  Subroutine to perform the decomposition of a vector, dvec = dmax * dmin,
     *  to ensure all elements greater than one are in dmax,
     *  and all elements less than one are in dmin.
     *  Input: dvec
     *  Output: dmax, dmin
     */

    assert(dvec.size() == dmax.size());
    assert(dvec.size() == dmin.size());

    const int ndim = dvec.size();
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

void mult_v_invd_u(const matXd& vmat, const vecXd& dvec, const matXd& umat, matXd& zmat) {
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

    const int ndim = vmat.rows();

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

void mult_v_d_u(const matXd& vmat, const vecXd& dvec, const matXd& umat, matXd& zmat) {
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

    const int ndim = vmat.rows();

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

void compute_Green(const matXd& U, const vecXd& S, const matXd& V, matXd& gtt, matXd& gt0, bool bool_measure_dynamic) {
    /*
     * returns (1 + USV^T)^-1, with method of QR decomposition
     * if bool_measure_dynamic, return (1 + USV^T)^-1 * USV^T too.
     */

    // split S = Sbi^-1 * Ss
    vecXd Sbi(S.size());
    vecXd Ss(S.size());
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
    matXd H = Sbi.asDiagonal() * U.transpose() + Ss.asDiagonal() * V.transpose();

    /* gtt */
    gtt = H.fullPivHouseholderQr().solve(Sbi.asDiagonal()*U.transpose());

    /* gt0 */
    if (bool_measure_dynamic) {
        gt0 = H.fullPivHouseholderQr().solve(Ss.asDiagonal()*V.transpose());
    }
}

void compute_Green(const SvdStack& left, const SvdStack& right, matXd& gtt, matXd& gt0, bool bool_measure_dynamic) {
    /*
     *  returns (1 + left * right^T)^-1 in a stable manner, with method of MGS factorization
     *  note: (1 + left * right^T)^-1 = (1 + (USV^T)_left * (VSU^T)_right)^-1
     */
    assert(left.n == right.n);
    const int ndim = left.n;

    /* at l = 0 */
    if(left.empty()) {
        compute_Green(right.matrixV(), right.singularValues(), right.matrixU(), gtt, gt0, bool_measure_dynamic);
        return;
    }

    /* at l = lt */
    if(right.empty()) {
        compute_Green(left.matrixU(), left.singularValues(), left.matrixV(), gtt, gt0, bool_measure_dynamic);
        return;
    }

    // local params
    const matXd& ul = left.matrixU();
    const vecXd& dl = left.singularValues();
    const matXd& vl = left.matrixV();
    const matXd& ur = right.matrixU();
    const vecXd& dr = right.singularValues();
    const matXd& vr = right.matrixV();

    vecXd dlmax(dl.size()), dlmin(dl.size());
    vecXd drmax(dr.size()), drmin(dr.size());

    matXd Atmp(ndim, ndim), Btmp(ndim, ndim);

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

    matXd tmp = Atmp + Btmp;
    mult_v_invd_u(ur, drmax, tmp.inverse(), Atmp);

    /* gtt */
    mult_v_invd_u(Atmp, dlmax, ul.transpose(), gtt);

    /* gt0 */
    if (bool_measure_dynamic) {
        mult_v_d_u(Atmp, dlmin, vl.transpose(), gt0);
    }
}

#endif //DQMC_HUBBARD_STABLEGREENS_HPP