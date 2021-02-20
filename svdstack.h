#ifndef HUBBARD_V1_3_SVDSTACK_H
#define HUBBARD_V1_3_SVDSTACK_H

#include <iostream>
#include <cmath>
#include <cassert>
#include <random>

#include <Eigen/LU>
#include <Eigen/SVD>
#include <unsupported/Eigen/MatrixFunctions>

typedef Eigen::MatrixXd matXd;
typedef Eigen::VectorXd vecXd;


/** UDV stack of a matrix product: U * D * V^T = ... A_2 * A_1 * A_0 */
struct SvdStack {

    typedef Eigen::JacobiSVD<matXd, Eigen::NoQRPreconditioner> SVD;
    std::vector<SVD> stack;
    int n;
    matXd tmp, V;
    int len = 0;


    SvdStack() = default;
    SvdStack(int n, int l)
            : n(n), tmp(n, n), V(n, n)
    {
        stack.reserve(l);
        for(int i = 0; i < l; ++i)
            stack.emplace_back(n, n, Eigen::ComputeFullU | Eigen::ComputeFullV);
    }

    bool empty() const
    {
        return len == 0;
    }

    /** prepends a matrix to the decomposition */
    void push(const matXd& m)
    {
        assert( m.rows() == n && m.cols() == n);
        assert( len < stack.size());

        if (len == 0) {
            stack[len].compute(m);
            V = stack[len].matrixV();
        }
        else {
            /** IMPORTANT! Mind the order of multiplication!
             *  Avoid confusing of different eigen scales here */
            tmp = ( m * matrixU() ) * singularValues().asDiagonal();
            stack[len].compute(tmp);
            V = V * stack[len].matrixV();
        }

        len += 1;
    }

    const vecXd& singularValues() const
    {
        assert(len > 0);
        return stack[len-1].singularValues();
    }

    const matXd& matrixU() const
    {
        assert(len > 0);
        return stack[len-1].matrixU();
    }

    matXd matrixV() const
    {
        assert(len > 0);
        matXd r = stack[0].matrixV();
        for (int i = 1; i < len; ++i) {
            r = r * stack[i].matrixV();
        }
        return r;
    }

    void pop()
    {
        assert(len > 0);
        len -= 1;
        if (len > 0) { V = V * stack[len].matrixV().transpose();}
    }

    void clear()
    {
        len = 0;
    }
};

#endif //HUBBARD_V1_3_SVDSTACK_H
