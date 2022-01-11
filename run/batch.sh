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
cb="false"
warm_up="true"
equal_measure="true"
dynamic_measure="true"

# generate folder name according to input params
if [ `echo "$u < 0.0"|bc` -eq 1 ] ; then
    u=$(echo "0.0 - $u"|bc)
    out_folder="L"$ll"b"${beta%.*}"u"${u%.*}
    u=$(echo "0.0 - $u"|bc)
else
    out_folder="L"$ll"b"${beta%.*}"u"${u%.*}
fi
# out_folder="example"

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
nwrap=$nwrap,nbin=$nbin,nsweep=$nsweep,cb=$cb,\
warm_up=$warm_up,equal_measure=$equal_measure,dynamic_measure=$dynamic_measure,\
out_folder_path=$out_folder_path ./run.sh

exit 0