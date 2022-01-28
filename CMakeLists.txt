cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 11)
project(dqmc_hubbard)

set(CMAKE_CXX_FLAGS "-march=native -O3 -fopenmp")

set(MKLROOT /public1/soft/oneAPI/2022.1/mkl/latest)
set(MPIROOT /public1/soft/oneAPI/2022.1/mpi/latest)
set(CMAKE_CXX_COMPILER ${MPIROOT}/bin/mpicxx)
set(CMAKE_C_COMPILER ${MPIROOT}/bin/mpicc)

set(BOOSTROOT /public1/soft/boost/1-55-0)
set(Boost_USE_RELEASE_lIBS ON)
set(Boost_USE_MULTITHREAD ON)

set(EIGENROOT /public1/home/sc81233/xuchengwang/local/eigen/eigen-3.4.0)

include_directories(
        ${MKLROOT}/include
        ${MPIROOT}/include
        ${BOOSTROOT}/include/boost
        ${EIGENROOT}
)

link_directories(
        ${MKLROOT}/lib/intel64
        ${MPIROOT}/lib
        ${BOOSTROOT}/lib
)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC_DIR)

add_executable(${PROJECT_NAME} ${SRC_DIR})

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)

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
        boost_mpi
        boost_serialization
)
