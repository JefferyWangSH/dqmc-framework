#!/bin/bash

# load environment variables
module load cmake/3.21.2 

cmake -DCMAKE_BUILD_TYPE=Release -G "CodeBlocks - Unix Makefiles" ..
# cmake -DCMAKE_BUILD_TYPE=Debug -G "CodeBlocks - Unix Makefiles" ..

# # clear up environment variables
# module unload cmake/3.21.2

exit 0