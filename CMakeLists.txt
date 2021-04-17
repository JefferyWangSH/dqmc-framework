cmake_minimum_required(VERSION 3.16)

project(dqmc_hubbard)

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -march=native -O3 -fopenmp" )

set(SOURCES
    src/hubbard.h
    src/hubbard.cpp
    src/detQMC.h
    src/detQMC.cpp
    src/svdstack.hpp
    src/StableGreens.hpp
    src/ProgressBar.hpp
)

# link eigen to MKL
include_directories(/opt/intel/compilers_and_libraries_2018/linux/mkl/include)

include_directories(/opt/intel/mkl/include)

include_directories(/usr/include/eigen3)

include_directories(${EIGEN3_INCLUDE_DIR})

link_directories(/opt/intel/mkl/lib/intel64)

link_directories(/opt/intel/lib/intel64)

add_executable(${PROJECT_NAME} src/main.cpp ${SOURCES})

target_link_libraries(${PROJECT_NAME} mkl_gf_lp64 mkl_intel_thread mkl_gnu_thread mkl_blas95_lp64 mkl_intel_lp64 mkl_sequential iomp5 mkl_core pthread m)

# boost
set(Boost_USE_RELEASE_lIBS ON)

set(Boost_USE_MULTITHREAD ON)

include_directories(/usr/local/include/boost)

link_directories(/usr/local/lib)

target_link_libraries(${PROJECT_NAME} boost_program_options)