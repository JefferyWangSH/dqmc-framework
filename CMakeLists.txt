cmake_minimum_required(VERSION 3.20)

project(dqmc_hubbard)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -march=native -O3 -fopenmp")

## set up boost
set(Boost_USE_RELEASE_lIBS ON)
set(Boost_USE_MULTITHREAD ON)

include_directories(
        ${EIGEN3_INCLUDE_DIR}
        /usr/include/eigen3
        /opt/intel/oneapi/mkl/latest/include/
        /usr/local/include/boost
)

link_directories(
        /opt/intel/oneapi/mkl/latest/lib/intel64                        # mkl
        /opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64    # linux openmp library, specifically libiomp5.so
        /usr/local/lib                                                  # boost
)

add_executable(${PROJECT_NAME}
        src/main.cpp
        src/Hubbard.cpp
        src/DetQMC.cpp
        src/SvdStack.cpp
        src/CheckerBoard.cpp
        src/BMatrixMult.hpp
        src/StableGreens.hpp
        src/ProgressBar.hpp
        src/MeasureData.cpp
        src/EqtimeMeasure.cpp
        src/EqtimeMeasureLinear.cpp
        src/DynamicMeasure.cpp
        src/DynamicMeasureLinear.cpp
)

target_include_directories(${PROJECT_NAME}
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(${PROJECT_NAME}
        mkl_gf_lp64
        mkl_intel_thread
        mkl_gnu_thread
        mkl_blas95_lp64
        mkl_intel_lp64
        mkl_sequential
        iomp5
        mkl_core
        pthread
        m
        boost_program_options
)
