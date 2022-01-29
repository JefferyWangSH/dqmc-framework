#ifndef DQMC_HUBBARD_SVDSTACK_H
#define DQMC_HUBBARD_SVDSTACK_H

/**
  *  This head file includes svd and SvdStack class for stable multiplication of long chains of matrices.
  *  BLAS and LAPACK libraries are needed for the svd decomposition.
  */

#include <vector>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


/** Structure for svd  */
class Svd {
private:
    Eigen::MatrixXd _u;
    Eigen::VectorXd _s;
    Eigen::MatrixXd _v;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Svd() = default;
    explicit Svd(int dim): _u(dim, dim), _s(dim), _v(dim, dim) {}

    Eigen::MatrixXd& MatrixU() { return _u; }
    Eigen::VectorXd& SingularValues() { return _s; }
    Eigen::MatrixXd& MatrixV() { return _v; }
};


/** UDV stack of a matrix product: U * D * V^T = ... A_2 * A_1 * A_0 */
class SvdStack {

public:
    std::vector<Svd> _stack;
    int _dim{1};
    int _len{0};
    Eigen::MatrixXd _tmp_matrix{};

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /* Construct functions */
    SvdStack() = default;

    SvdStack(int dim, int len);

    bool empty() const;

    int dim() const;

    int length() const;

    void clear();

    /** prepends a matrix to the decomposition */
    void push(const Eigen::MatrixXd &matrix);

    /** pop the last matrix from the decomposition */
    void pop();

    const Eigen::VectorXd SingularValues();

    const Eigen::MatrixXd MatrixU();

    const Eigen::MatrixXd MatrixV();

};

#endif //DQMC_HUBBARD_SVDSTACK_H