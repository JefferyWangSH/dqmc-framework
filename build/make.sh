#!/bin/bash

# load environment variables
module purge
module load gcc/10.2.0
module load oneAPI/2022.1
module load mpi/intel/2022.1

# build
make

module purge
exit 0
