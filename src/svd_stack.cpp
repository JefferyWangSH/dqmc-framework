#include <iostream>
#include <cmath>
#include <cassert>
#include "svd_stack.h"
#include "mkl_lapacke.h"
#include "utils/linear_algebra.hpp"


namespace Utils {

    using VecSvd = std::vector<SvdClass>;
    using Matrix = Eigen::MatrixXd;
    using Vector = Eigen::VectorXd;

    SvdStack::SvdStack(int mat_dim, int stack_length) 
                : m_mat_dim(mat_dim), 
                  m_tmp_matrix(mat_dim, mat_dim)
    {
        this->m_stack.reserve(stack_length);
        for (int i = 0; i < stack_length; ++i) {
            this->m_stack.emplace_back(mat_dim);
        }
    }

    bool SvdStack::empty() const { return this->m_stack_length == 0; }

    int SvdStack::MatDim() const { return this->m_mat_dim; }

    int SvdStack::StackLength() const { return this->m_stack_length; }

    void SvdStack::clear() { this->m_stack_length = 0; }


    void SvdStack::push(const Matrix& matrix) {
        assert( matrix.rows() == this->m_mat_dim && matrix.cols() == this->m_mat_dim );
        assert( this->m_stack_length < (int)this->m_stack.size() );

        if (this->m_stack_length == 0) {
            // udv decomposition
            Utils::LinearAlgebra::mkl_lapack_dgesvd (
                this->m_mat_dim, 
                this->m_mat_dim, 
                matrix, 
                this->m_stack[this->m_stack_length].MatrixU(), 
                this->m_stack[this->m_stack_length].SingularValues(), 
                this->m_stack[this->m_stack_length].MatrixV() );
        }
        else {
            // important! mind the order of multiplication!
            // Avoid mixing of different numerical scales here
            this->m_tmp_matrix = ( matrix * this->MatrixU() ) * this->SingularValues().asDiagonal();
            Utils::LinearAlgebra::mkl_lapack_dgesvd (
                this->m_mat_dim, 
                this->m_mat_dim, 
                this->m_tmp_matrix, 
                this->m_stack[this->m_stack_length].MatrixU(), 
                this->m_stack[this->m_stack_length].SingularValues(), 
                this->m_stack[this->m_stack_length].MatrixV() );
        }
        this->m_stack_length += 1;
    }

    void SvdStack::pop() {
        // caution that the memory is not actually released
        assert(this->m_stack_length > 0);
        this->m_stack_length -= 1;
    }


    const Vector SvdStack::SingularValues() {
        assert(this->m_stack_length > 0);
        return this->m_stack[this->m_stack_length-1].SingularValues();
    }

    const Matrix SvdStack::MatrixU() {
        assert(this->m_stack_length > 0);
        return this->m_stack[this->m_stack_length-1].MatrixU();
    }

    const Matrix SvdStack::MatrixV() {
        assert(this->m_stack_length > 0);
        Matrix r = this->m_stack[0].MatrixV();
        for (int i = 1; i < this->m_stack_length; ++i) {
            r = r * this->m_stack[i].MatrixV();
        }
        return r;
    }


} // namespace Utils
