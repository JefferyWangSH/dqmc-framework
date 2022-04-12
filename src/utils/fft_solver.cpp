#include "utils/fft_solver.h"

namespace Utils {

    namespace FFTSolver {

        void FFTSolver2d::set_up_dimension(int row, int col) {
            this->m_row = row;
            this->m_col = col;
        }

        void FFTSolver2d::initial() {
            // allocate for intermediate variable
            if (this->m_dft_data) { this->m_dft_data.reset(); }
            this->m_dft_data = std::make_unique<std::complex<double>[]>(this->m_row*this->m_col);
                    
            // create and commit fft descriptor                
            MKL_LONG dim_sizes[2] = {this->m_row, this->m_col};
            DftiCreateDescriptor(&this->m_desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, dim_sizes);
            DftiCommitDescriptor(this->m_desc_handle);
        }

        void FFTSolver2d::compute(const Eigen::MatrixXd &in, Eigen::MatrixXd &out) {
            // TODO: accelerate by using r2c fft
            // the input matrix should be real
            assert( in.rows() == this->m_row );
            assert( in.cols() == this->m_col );

            // map input eigen matrix to c-style array
            Eigen::Map<Eigen::MatrixXcd>(&this->m_dft_data[0], this->m_row, this->m_col) = in;

            // in-place fft
            DftiComputeForward(this->m_desc_handle, &this->m_dft_data[0]);

            // map transformed data to eigen matrix
            // output only the eal part of the complex results 
            out = Eigen::Map<Eigen::MatrixXcd>(&this->m_dft_data[0], this->m_row, this->m_col).real();
        }

        void FFTSolver2d::deallocate() {
            // deallocate memory
            if (this->m_dft_data) { this->m_dft_data.reset(); }
            DftiFreeDescriptor(&m_desc_handle);
        }

    } // namespace FFTSolver

} // namespace Utils