#ifndef DQMC_HUBBARD_FFT_SOLVER_H
#define DQMC_HUBBARD_FFT_SOLVER_H
#pragma once

// This class, which is initially intended to transform real-space greens function to momentum space, 
// is temporarily not used in the dqmc simulations.
// This is because that the transformation in dqmc is actually a four-dimension dft projected to 
// a two-dimensional surface of momentum (kx, ky) due to the conservation of momentum
//
//     < c(k1) c^\dagger(k2) >    ->    k1x = k2x and k1y = k2y 
//
// for ki = (kix, kiy) being two-dimensional momentum.
// which should be more straightforward and efficient when intuitively calculated by `for` loops.

/**
  *  This header file includes `FFTSolver2d` class for the discrete fourier
  *  transformation of real-space data using intel mkl fft, 
  *  with friendly interface of Eigen::Matrix.
  */

#include <complex>
#include <memory>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include <mkl_dfti.h>


namespace FFTSolver {

    /**
      *  Specialized fft solver class for two-dimensional fft of real input data
      *  with interface of Eigen::MatrixXd class.
      *  Note: In principle, the output data should be in complex form, while in the usage of FFTSolver 
      *        only the real part of the transformed results would be conserved (r2c). 
      *        This is mainly because only the real part of transformation is of physical interpretations in dqmc simulations. 
      *  Todo: accelarate the transformation using r2cfft directly ( considering the symmetry )
      */
    class FFTSolver2d {
        private:
            // dimension information
            int row{}, col{};

            // variables storing the intermediate data during fft
            std::unique_ptr<std::complex<double>[]> dft_data{};

            // descriptor of mkl fft
            DFTI_DESCRIPTOR_HANDLE desc_handle{NULL};

        public:
            // default constructing function
            FFTSolver2d() = default;

            // set up dimensions of fft
            void set_up_dimension(int row, int col);

            // initialization
            void initial();

            // fft computation
            void compute(const Eigen::MatrixXd &in, Eigen::MatrixXd &out);

            // deallocate memory
            void deallocate();
    };

} // namespace FFTSolver

#endif  // DQMC_HUBBARD_FFT_SOLVER_H
