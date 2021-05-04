#include <iostream>
#include <cmath>
#include <cassert>
#include "SvdStack.h"


SvdStack::SvdStack(int n, int l): n(n), tmp(n, n) {
    stack.reserve(l);
    for(int i = 0; i < l; ++i) {
        stack.emplace_back(n);
    }
}

void SvdStack::resize(int n_, int l_) {
    SvdStack newStack(n_, l_);
    *this = newStack;
}

bool SvdStack::empty() const {
    return len == 0;
}

void SvdStack::clear() {
    len = 0;
}


/** svd decomposition of arbitrary M * N real matrix, using LAPACK */
static void lapack_dgesdd(const int &m, const int &n, const Eigen::MatrixXd &a, Eigen::MatrixXd &u, Eigen::VectorXd &s, Eigen::MatrixXd &v) {
    assert(m == a.rows());
    assert(n == a.cols());

    /* Matrix size */
    const int lda = m;
    const int ldu = m;
    const int ldvt = n;

    int info, lwork;
    double wkopt;
    double* work;

    /* Local arrays */
    /* iwork dimension should be at least 8*min(m,n) */
    int iwork[8 * n];
    double s_[n], u_[ldu * m], vt_[ldvt * n];
    double a_[lda * n];

    for (int i = 0; i < lda * n; ++i) {
        a_[i] = a(i / lda, i % lda);
    }

    /* Query and allocate the optimal workspace */
    lwork = -1;
    dgesdd( "A", &m, &n, a_, &lda, s_, u_, &ldu, vt_, &ldvt, &wkopt, &lwork, iwork, &info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );

    /* Compute SVD */
    dgesdd( "A", &m, &n, a_, &lda, s_, u_, &ldu, vt_, &ldvt, work, &lwork, iwork, &info );

    /* Check for convergence */
    if( info > 0 ) {
        std::cerr << "The algorithm computing SVD failed to converge." << std::endl;
        exit( 1 );
    }

    /* Convert to Eigen */
    u = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(vt_, n, n);
    s = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(s_, 1, lda);
    v = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(u_, m, m);

    /* Free workspace */
    free( (void*)work );
}


void SvdStack::push(const Eigen::MatrixXd &a) {
    assert( a.rows() == n && a.cols() == n );
    assert( len < stack.size() );

    if (len == 0) {
        lapack_dgesdd(n, n, a, stack[len].matrixU(), stack[len].singularValues(), stack[len].matrixV());
    }
    else {
        /** IMPORTANT! Mind the order of multiplication!
         *  Avoid confusing of different eigen scales here */
        tmp = ( a * matrixU() ) * singularValues().asDiagonal();
        lapack_dgesdd(n, n, tmp, stack[len].matrixU(), stack[len].singularValues(), stack[len].matrixV());
    }
    len += 1;
}

void SvdStack::pop() {
    assert(len > 0);
    len -= 1;
}

Eigen::VectorXd SvdStack::singularValues() {
    assert(len > 0);
    return stack[len-1].singularValues();
}

Eigen::MatrixXd SvdStack::matrixU() {
    assert(len > 0);
    return stack[len-1].matrixU();
}

Eigen::MatrixXd SvdStack::matrixV(){
    assert(len > 0);
    Eigen::MatrixXd r = stack[0].matrixV();
    for (int i = 1; i < len; ++i) {
        r = r * stack[i].matrixV();
    }
    return r;
}
