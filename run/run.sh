#!/bin/bash

#SBATCH --partition=v6_384
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --job-name=dqmc
#SBATCH --output=../results/log.log

exe="../build/dqmc_hubbard"
ll=4
lt=80
b=4.0
u=-4.0
mu=0.0
nbin=24
nsweep=100
cb="false"
equal_measure="true"
dynamic_measure="true"

mpirun ${exe} --ll=${ll} --lt=${lt} --beta=${b} --u=${u} --mu=${mu} --checkerboard=${cb} --eqtime=${equal_measure} --dynamic=${dynamic_measure} --nbin=${nbin} --nsweep=${nsweep}

exit 0