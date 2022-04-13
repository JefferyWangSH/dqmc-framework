#ifndef UTIL_SVD_STACK_H
#define UTIL_SVD_STACK_H

/**
  *  This head file includes SvdClass and SvdStack class for stable 
  *  multiplication of long chains of dense matrices.
  *  BLAS and LAPACK libraries are needed for the svd decomposition.
  */

#include <vector>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace Utils {


    // --------------------------------------- Utils::SvdClass ------------------------------------------
    class SvdClass {
        private:
            using uMat = Eigen::MatrixXd;
            using sVec = Eigen::VectorXd;
            using vMat = Eigen::MatrixXd;
            
            uMat m_u_mat{};
            sVec m_s_vec{};
            vMat m_v_mat{};

        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            SvdClass() = default;

            explicit SvdClass(int dim): m_u_mat(dim, dim), m_s_vec(dim), m_v_mat(dim, dim) {}

            uMat& MatrixU() { return this->m_u_mat; }
            sVec& SingularValues() { return this->m_s_vec; }
            vMat& MatrixV() { return this->m_v_mat; }
    };


    // ---------------- Utils::SvdStack class for multiplications of matrix stacks ----------------------
    // udv stack of a matrix product: u * d * vt = ... A_2 * A_1 * A_0
    class SvdStack {
        private:
            using VecSvd = std::vector<SvdClass>;
            using Matrix = Eigen::MatrixXd;
            using Vector = Eigen::VectorXd;
        
        public:
            VecSvd m_stack{};
            int m_mat_dim{};  
            int m_stack_length{0};

            Matrix m_tmp_matrix{};

        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            SvdStack() = default;

            explicit SvdStack(int mat_dim, int stack_length);

            // interface
            bool empty() const;
            int MatDim() const;
            int StackLength() const;

            // return udv decomposition matrices of the stack
            const Vector SingularValues();
            const Matrix MatrixU();
            const Matrix MatrixV();
            
            // clear the stack
            // simply set stack_length = 0, the memory is not deallocated.
            void clear();

            // defined operations of the class: push and pop
            // adding a matrix to the stack from the left
            void push(const Matrix& matrix);

            // pop the last matrix from the stack
            void pop();

    };

} // namespace Utils

#endif // UTIL_SVD_STACK_H