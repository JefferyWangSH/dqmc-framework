#!/bin/bash

# load environment variables
module purge
module load cmake/3.21.2 
module load gcc/10.2.0
module load oneAPI/2022.1
module load mpi/intel/2022.1

cmake -DCMAKE_BUILD_TYPE=Release -G "CodeBlocks - Unix Makefiles" ..
# cmake -DCMAKE_BUILD_TYPE=Debug -G "CodeBlocks - Unix Makefiles" ..

# clear up environment variables
module purge
exit 0