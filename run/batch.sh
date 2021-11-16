#!/bin/bash

# name of output folder
folder_name="example"

# create output folder if not exist
folder_path="../results/"${folder_name}
if [ ! -d ${folder_path} ]; then
  mkdir ${folder_path}
fi

# set up jobname and log output name
jobname=${folder_name}
output=${folder_path}"/log.log"

sbatch --job-name ${jobname} --output ${output} ./run.sh

exit 0