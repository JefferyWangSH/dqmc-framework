#!/bin/bash

# parmas for submitting
partition="v6_384"
nodes=1
ntasks_per_node=20
cpus_per_task=1

# customized params
exe="../build/dqmc_hubbard"
ll=4
lt=160
beta=8.0
u=-4.0
mu=-0.2
nwrap=10
nbin=20
nsweep=100
checker_board="false"
warm_up="true"
equal_measure="true"
dynamic_measure="true"

# Supported physical observables:
#
#   1. filling_number                   (equal-time)
#   2. double_occupancy                 (equal-time)
#   3. kinetic_energy                   (equal-time)
#   4. momentum_distribution            (equal-time)
#   5. local_spin_corr                  (equal-time)
#   6. spin_density_structure_factor    (equal-time)
#   7. charge_density_structure_factor  (equal-time)
#   8. matsubara_greens                 (dynamic)
#   9. density_of_states                (dynamic)
#  10. superfluid_stiffness             (dynamic)
#
# options should be separated by space
# option 'all' : measure all supported observables
# option 'none': no measurement
obs_list="all"

# # generate folder name according to input params
# if [ `echo "$u < 0.0"|bc` -eq 1 ] ; then
#     u=$(echo "0.0 - $u"|bc)
#     out_folder="L"$ll"b"${beta%.*}"u"${u%.*}
#     u=$(echo "0.0 - $u"|bc)
# else
#     out_folder="L"$ll"b"${beta%.*}"u"${u%.*}
# fi
out_folder="example"

# create output folder if not exist
out_folder_path="../results/"${out_folder}
if [ ! -d ${out_folder_path} ]; then
  mkdir -p ${out_folder_path}
fi

# set up jobname and log output name
jobname=$out_folder
output=$out_folder_path"/log.log"
error=$out_folder_path"/err.log"

# set up new environment and pass variables to scripts
sbatch --job-name=${jobname} --output=${output} --error=${error} \
--partition=$partition --nodes=$nodes --ntasks-per-node=$ntasks_per_node --cpus-per-task=$cpus_per_task \
--export=exe=$exe,ll=$ll,lt=$lt,beta=$beta,u=$u,mu=$mu,\
nwrap=$nwrap,nbin=$nbin,nsweep=$nsweep,checker_board=$checker_board,\
warm_up=$warm_up,equal_measure=$equal_measure,dynamic_measure=$dynamic_measure,\
obs_list="$obs_list",out_folder_path=$out_folder_path ./run.sh

exit 0