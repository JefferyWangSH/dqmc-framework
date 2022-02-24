#!/bin/bash

# load environment variables
module purge
module load cmake/3.21.2 
module load gcc/10.2.0
module load oneAPI/2022.1
module load mpi/intel/2022.1

cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -G "CodeBlocks - Unix Makefiles" ..
# cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -G "CodeBlocks - Unix Makefiles" ..

# clear up environment variables
module purge
exit 0