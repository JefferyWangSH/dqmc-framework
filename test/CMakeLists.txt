# project
cmake_minimum_required(VERSION 3.21)
project(test)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-march=native -O3 -fopenmp -Wall")
get_filename_component(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# add targets and links
file(GLOB SOURCE_FILE 
    # ${PROJECT_SOURCE_DIR}/*/*.cpp ${PROJECT_SOURCE_DIR}/*/*.hpp 
    ${PROJECT_SOURCE_DIR}/test/test_main.cpp
    ${PROJECT_SOURCE_DIR}/*/*/*.cpp  )
# list(REMOVE_ITEM SOURCE_FILE ${PROJECT_SOURCE_DIR}/src/dqmc_main.cpp)
add_executable(
    ${PROJECT_NAME}  ${SOURCE_FILE} 
    ${PROJECT_SOURCE_DIR}/src/dqmc.cpp
    ${PROJECT_SOURCE_DIR}/src/dqmc_walker.cpp
    ${PROJECT_SOURCE_DIR}/src/svd_stack.cpp
    ${PROJECT_SOURCE_DIR}/src/fft_solver.cpp
    ${PROJECT_SOURCE_DIR}/src/random.cpp
    )
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
find_package(Eigen3 MODULE 3.3.7 REQUIRED)
if (Eigen3_FOUND)
    message(STATUS "Found Eigen3 (eigen3_include_dir): ${EIGEN3_INCLUDE_DIR}")
    target_link_libraries (${PROJECT_NAME} PRIVATE Eigen3::Eigen)
else()
    message(FATAL_ERROR "Eigen3 not found")
endif()

# find Boost
find_package(Boost MODULE 1.71.0 COMPONENTS program_options mpi serialization REQUIRED)
if (Boost_FOUND)
    message(STATUS "Found Boost: version ${Boost_VERSION}")
    message(STATUS "Found Boost (boost_include_dirs): ${Boost_INCLUDE_DIRS}")
    message(STATUS "Found Boost (boost_library_dirs): ${Boost_LIBRARY_DIRS}")
    set(Boost_USE_RELEASE_lIBS ON)
    set(Boost_USE_MULTITHREAD ON)
    target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIRS})
    target_link_libraries (${PROJECT_NAME} PRIVATE ${Boost_LIBRARIES})
else()
    message(FATAL_ERROR "Boost not found")
endif()