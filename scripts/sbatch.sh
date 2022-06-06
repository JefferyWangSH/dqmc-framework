#!/bin/bash

########################################################### submitting params
partition="v6_384"
nodes=1
ntasks_per_node=10
cpus_per_task=1
mem_per_cpu=16G


########################################################### program params
# executable object
exe="../build/dqmc"

# program options
# edit the config file
config_file="../example/config.toml"

# model params
# for Repulsive and Attractive Hubbard model
hopping_t=1.0
onsite_u=4.0
chemical_potential=0.0
sed -i "s/hopping_t = .*/hopping_t = "$hopping_t"/" $config_file
sed -i "s/onsite_u = .*/onsite_u = "$onsite_u"/" $config_file
sed -i "s/chemical_potential = .*/chemical_potential = "$chemical_potential"/" $config_file

# lattice params
# side length of the lattice ( square and cubic )
l=4
# # for square lattice
# sed -i "s/cell = .*/cell = [ "$l", "$l" ]/" $config_file
# for cubic lattice
sed -i "s/cell = .*/cell = [ "$l", "$l", "$l" ]/" $config_file

# momentum for measurments of observables
momentum="RPoint"
momentum_list="KstarsAll"
sed -i "s/momentum = .*/momentum = \""$momentum"\"/" $config_file
sed -i "s/momentum_list = .*/momentum_list = \""$momentum_list"\"/" $config_file

# monte carlo params
# inverse temperature, imaginary-time grids and pace of stabilization
beta=8.0
time_size=160
stabilization_pace=10
sed -i "s/beta = .*/beta = "$beta"/" $config_file
sed -i "s/time_size = .*/time_size = "$time_size"/" $config_file
sed -i "s/stabilization_pace = .*/stabilization_pace = "$stabilization_pace"/" $config_file

# measuring params
# MC sweeps for warmup, number of bins, samples in one bin,
# and MC sweeps between adjacent bins 
sweeps_warmup=512
bin_num=20
bin_size=100
sweeps_between_bins=20
sed -i "s/sweeps_warmup = .*/sweeps_warmup = "$sweeps_warmup"/" $config_file
sed -i "s/bin_num = .*/bin_num = "$bin_num"/" $config_file
sed -i "s/bin_size = .*/bin_size = "$bin_size"/" $config_file
sed -i "s/sweeps_between_bins = .*/sweeps_between_bins = "$sweeps_between_bins"/" $config_file

# edit observable list
# sed -i "s/observables = .*/observables = [ \"filling_number\", \"double_occupancy\" ]/" $config_file
sed -i "s/observables = .*/observables = [ \"all\" ]/" $config_file

# set up name of the output folder 
folder_name="L"$l"U"$onsite_u"b"$beta
output_folder="../example/"$folder_name
fields_file=$output_folder"/fields.out"

# set up jobname and log output name
job_name=$folder_name
log_file=$output_folder"/log.out"
err_file=$output_folder"/err.out"

# create output folder if not exist
if [ ! -d $output_folder ]; then
    mkdir -p $output_folder
fi
# delete previous log file if exist
if [ -f $log_file ]; then
    rm $log_file && touch $log_file
fi


########################################################### submit mission
sbatch \
--job-name=$job_name --output=$log_file --error=$err_file --partition=$partition \
--nodes=$nodes --ntasks-per-node=$ntasks_per_node --cpus-per-task=$cpus_per_task --mem-per-cpu=$mem_per_cpu \
--export=exe=$exe,config_file=$config_file,output_folder=$output_folder,fields_file=$fields_file \
./run.sh


exit 0