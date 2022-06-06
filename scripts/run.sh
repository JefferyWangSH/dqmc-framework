#!/bin/bash

# load modules
module purge 
module load gcc/10.2.0
module load oneAPI/2022.1
module load mpi/intel/2022.1

if [ -f $fields_file ]; then
    # input fields configs detected
    mpirun $exe --config $config_file --fields $fields_file --output $output_folder

else 
    # simulating with random fields configs
    mpirun $exe --config $config_file --output $output_folder

fi

module purge

exit 0