```bash
module use ~/modules
module load vmd/1.9.4lm
conda activate vmdplugin # to ensure conda prefix in the right place 
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
vmd
```
