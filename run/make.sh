#!/bin/bash

module purge
# module load gcc/10.2.0
module load python/3.9.6
module load oneAPI/2022.1
module load mpi/intel/2022.1

cd ../build && make

module purge
exit 0