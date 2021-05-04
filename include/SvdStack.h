#ifndef DQMC_HUBBARD_SVDSTACK_H
#define DQMC_HUBBARD_SVDSTACK_H

/*
 *  This head file includes svd and SvdStack class for stable multiplication of long chains of matrices.
 *  BLAS and LAPACK libraries are needed for the svd decomposition.
 */

#include <vector>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


/** Structure for svd  */
class svd {
private:
    Eigen::MatrixXd u;
    Eigen::VectorXd s;
    Eigen::MatrixXd v;

public:
    svd() = default;
    explicit svd(int n): u(), s(n), v(n, n) {}

    Eigen::MatrixXd& matrixU() { return u; }
    Eigen::VectorXd& singularValues() { return s; }
    Eigen::MatrixXd& matrixV() { return v; }
};


/** UDV stack of a matrix product: U * D * V^T = ... A_2 * A_1 * A_0 */
class SvdStack
{
public:

    std::vector<svd> stack;
    int n{};
    Eigen::MatrixXd tmp{};
    int len = 0;


    /* Construct functions */
    SvdStack() = default;
    SvdStack(int n, int l);

    void resize(int n, int l);

    bool empty() const;

    void clear();

    /** prepends a matrix to the decomposition */
    void push(const Eigen::MatrixXd &mat);

    /** pop the last matrix from the decomposition */
    void pop();

    Eigen::VectorXd singularValues();

    Eigen::MatrixXd matrixU();

    Eigen::MatrixXd matrixV();

};

#endif //DQMC_HUBBARD_SVDSTACK_H