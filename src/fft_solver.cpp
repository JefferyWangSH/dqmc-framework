#include "fft_solver.h"

namespace FFTSolver {

    void FFTSolver2d::set_up_dimension(int row, int col) {
        this->row = row;
        this->col = col;
    }

    void FFTSolver2d::initial() {
        // allocate for intermediate variable
        if (this->dft_data) { this->dft_data.reset(); }
        this->dft_data = std::make_unique<std::complex<double>[]>(this->row*this->col);
                
        // create and commit fft descriptor                
        MKL_LONG dim_sizes[2] = {this->row, this->col};
        DftiCreateDescriptor(&this->desc_handle, DFTI_DOUBLE, DFTI_COMPLEX, 2, dim_sizes);
        DftiCommitDescriptor(this->desc_handle);
    }

    void FFTSolver2d::compute(const Eigen::MatrixXd &in, Eigen::MatrixXd &out) {
        // TODO: accelerate by using r2c fft
        // the input matrix should be real
        assert( in.rows() == this->row );
        assert( in.cols() == this->col );

        // map input eigen matrix to c-style array
        Eigen::Map<Eigen::MatrixXcd>(&this->dft_data[0], this->row, this->col) = in;

        // in-place fft
        DftiComputeForward(desc_handle, &this->dft_data[0]);

        // map transformed data to eigen matrix
        // output only the eal part of the complex results 
        out = Eigen::Map<Eigen::MatrixXcd>(&this->dft_data[0], this->row, this->col).real();
    }

    void FFTSolver2d::deallocate() {
        // deallocate memory
        if (this->dft_data) { this->dft_data.reset(); }
        DftiFreeDescriptor(&desc_handle);
    }

} // namespace FFTSolver