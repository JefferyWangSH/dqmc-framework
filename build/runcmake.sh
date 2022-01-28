#!/bin/bash

# load environment variables
module purge
module load cmake/3.21.2 
module load gcc/10.2.0

cmake -DCMAKE_BUILD_TYPE=Release -G "CodeBlocks - Unix Makefiles" ..

# clear up environment variables
module purge
exit 0