#!/bin/bash

#SBATCH --partition="v6_384"
#SBATCH --job-name="example"
#SBATCH --output="../example/log.out"
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=10

module load gcc/10.2.0
module load cmake/3.21.2
module load python/3.9.6
module load oneAPI/2022.1
module load mpi/intel/2022.1

mpirun ../build/dqmc --config ../example/config.toml --fields ../example/output/fields.out --output ../example/output

module purge
exit 0