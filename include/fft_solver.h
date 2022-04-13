#ifndef UTIL_FFT_SOLVER_H
#define UTIL_FFT_SOLVER_H
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


namespace Utils {

    namespace FFTSolver {

        // ---------------------------------- Utils::FFTSolver::FFTSolver2d ---------------------------------------
        /**
          *  Specialized fft solver class for two-dimensional fft of real input data
          *  with interface of Eigen::MatrixXd class.
          *  Note: In principle, the output data should be in complex form, while in the usage of FFTSolver 
          *        only the real part of the transformed results would be conserved (r2c). 
          *        This is mainly because only the real part of transformation is of physical 
          *        interpretations in dqmc simulations. 
          *  Todo: accelarate the transformation using r2cfft directly ( considering the symmetry )
          */
        class FFTSolver2d {
            private:
                using RowDim = int;
                using ColDim = int;
                using ptrCpxArray = std::unique_ptr<std::complex<double>[]>;
                using DftHandler = DFTI_DESCRIPTOR_HANDLE;
                using Matrix = Eigen::MatrixXd;


                // dimension information
                int m_row{}, m_col{};

                // variables storing the intermediate data during fft
                ptrCpxArray m_dft_data{};

                // descriptor of mkl fft
                DftHandler m_desc_handle{NULL};

            public:
                // default constructing function
                FFTSolver2d() = default;

                // set up dimensions of fft
                void set_up_dimension(int row, int col);

                // initialization
                void initial();

                // fft computation
                void compute(const Matrix& in, Matrix& out);

                // deallocate memory
                void deallocate();
        };

    } // namespace FFTSolver

} // namespace Utils


#endif // UTIL_FFT_SOLVER_H
