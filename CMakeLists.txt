# project
cmake_minimum_required(VERSION 3.21)
project(dqmc_hubbard)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-march=native -O3 -fopenmp")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# target
aux_source_directory(${PROJECT_SOURCE_DIR}/src SOURCE_DIR)
add_executable(${PROJECT_NAME} ${SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include)

# find MKL
find_package(MKL MODULE REQUIRED)
if (MKL_FOUND)
    message(STATUS "Found MKL (mkl_include_dir): ${MKL_INCLUDE_DIR}")
    message(STATUS "Found MKL (mkl_library_dir): ${MKL_LIBRARY_DIR}")
    target_include_directories(${PROJECT_NAME} PRIVATE ${MKL_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME} PRIVATE ${MKL_LIBRARIES})
else()
    message(FATAL_ERROR "MKL not found")
endif()

# find MPI
find_package(MPI MODULE REQUIRED)
if (MPI_CXX_FOUND)
    message(STATUS "Found MPI (mpi_cxx_include_path): ${MPI_CXX_INCLUDE_PATH}")
    target_include_directories(${PROJECT_NAME} PRIVATE ${MPI_CXX_INCLUDE_PATH})
    target_link_libraries (${PROJECT_NAME} PRIVATE MPI::MPI_CXX)
else()
    message(FATAL_ERROR "MPI not found")
endif()

# find Eigen3
find_package(Eigen3 MODULE 3.4.0 REQUIRED)
if (Eigen3_FOUND)
    message(STATUS "Found Eigen3 (eigen3_include_dir): ${EIGEN3_INCLUDE_DIR}")
    target_link_libraries (${PROJECT_NAME} PRIVATE Eigen3::Eigen)
else()
    message(FATAL_ERROR "Eigen3 not found")
endif()

# find Boost
find_package(Boost MODULE 1.55.0 COMPONENTS program_options mpi serialization REQUIRED)
if (Boost_FOUND)
    message(STATUS "Found Boost: version ${Boost_LIB_VERSION}")
    message(STATUS "Found Boost (boost_include_dirs): ${Boost_INCLUDE_DIRS}")
    message(STATUS "Found Boost (boost_libraries): ${Boost_LIBRARIES}")
    set(Boost_USE_RELEASE_lIBS ON)
    set(Boost_USE_MULTITHREAD ON)
    target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries (${PROJECT_NAME} PRIVATE ${Boost_LIBRARIES})
else()
    message(FATAL_ERROR "Boost not found")
endif()