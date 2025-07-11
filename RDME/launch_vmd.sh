#! /bin/bash

# This script is used to launch VMD

cd /projects/beyi/
module load 
source activate /projects/beyi/sw/conda/envs/vmdplugin
module load vmd/1.9.4lm
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
vmd /projects/beyi/$USER/SummerSchool_2025/RDME/trajectory/MinCell_1.lm -e /projects/beyi/$USER/SummerSchool_2025/RDME/trajectory/render.vmd