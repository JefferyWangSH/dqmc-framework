#!/bin/bash

mpirun $exe --ll=$ll --lt=$lt --beta=$beta --u=$u --mu=$mu \
--checkerboard=$cb --warm-up=$warm_up --eqtime-measure=$equal_measure --dynamic-measure=$dynamic_measure \
--nwrap=$nwrap --nbin=$nbin --nsweep=$nsweep --out-folder-path=$out_folder_path

exit 0