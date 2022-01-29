#include <iostream>
#include <cmath>
#include <cassert>
#include "mkl_lapacke.h"
#include "svd_stack.h"
#include "matrix_util.hpp"

SvdStack::SvdStack(int dim, int len): _dim(dim), _tmp_matrix(dim, dim) {
    this->_stack.reserve(len);
    for(int i = 0; i < len; ++i) {
        this->_stack.emplace_back(dim);
    }
}

bool SvdStack::empty() const {
    return this->_len == 0;
}

int SvdStack::dim() const {
    return this->_dim;
}

int SvdStack::length() const {
    return this->_len;
}

void SvdStack::clear() {
    this->_len = 0;
}

void SvdStack::push(const Eigen::MatrixXd &a) {
    assert( a.rows() == this->_dim && a.cols() == this->_dim );
    assert( this->_len < this->_stack.size() );

    if (this->_len == 0) {
        MatrixUtil::mkl_lapack_dgesvd(this->_dim, this->_dim, a, 
        this->_stack[this->_len].MatrixU(), this->_stack[this->_len].SingularValues(), this->_stack[this->_len].MatrixV());
    }
    else {
        /** IMPORTANT! Mind the order of multiplication!
         *  Avoid confusing of different eigen scales here */
        this->_tmp_matrix = ( a * this->MatrixU() ) * this->SingularValues().asDiagonal();
        MatrixUtil::mkl_lapack_dgesvd(this->_dim, this->_dim, this->_tmp_matrix, 
        this->_stack[this->_len].MatrixU(), this->_stack[this->_len].SingularValues(), this->_stack[this->_len].MatrixV());
    }
    this->_len += 1;
}

void SvdStack::pop() {
    assert(this->_len > 0);
    this->_len -= 1;
}

const Eigen::VectorXd SvdStack::SingularValues() {
    assert(this->_len > 0);
    return this->_stack[this->_len-1].SingularValues();
}

const Eigen::MatrixXd SvdStack::MatrixU() {
    assert(this->_len > 0);
    return this->_stack[this->_len-1].MatrixU();
}

const Eigen::MatrixXd SvdStack::MatrixV() {
    assert(this->_len > 0);
    Eigen::MatrixXd r = this->_stack[0].MatrixV();
    for (int i = 1; i < this->_len; ++i) {
        r = r * this->_stack[i].MatrixV();
    }
    return r;
}
