#include <iostream>
#include <cmath>
#include <cassert>

#include "mkl_lapacke.h"

#include "svd_stack.h"
#include "matrix_util.hpp"


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

void SvdStack::push(const Eigen::MatrixXd &a) {
    assert( a.rows() == n && a.cols() == n );
    assert( len < stack.size() );

    if (len == 0) {
        MatrixUtil::mkl_lapack_dgesvd(n, n, a, stack[len].matrixU(), stack[len].singularValues(), stack[len].matrixV());
    }
    else {
        /** IMPORTANT! Mind the order of multiplication!
         *  Avoid confusing of different eigen scales here */
        tmp = ( a * matrixU() ) * singularValues().asDiagonal();
        MatrixUtil::mkl_lapack_dgesvd(n, n, tmp, stack[len].matrixU(), stack[len].singularValues(), stack[len].matrixV());
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
