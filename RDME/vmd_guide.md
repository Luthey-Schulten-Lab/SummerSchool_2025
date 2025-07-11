# Visualizing 4DWCM Trajectory

This guide explains how to visualize 4DWCM model trajectories using [VMD 1.9.4a](https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD).

VMD and the [VMD plugin](https://github.com/Luthey-Schulten-Lab/LMVMDPlugin) are pre-installed on NCSA Delta HPC.

Since VMD requires a graphical interface, we'll use Open OnDemand's Desktop interactive app to access a Linux GUI for launching VMD and viewing trajectories.

## 1. Initialize the OOD Interactive Session
1. Navigate to the Open OnDemand dashboard.

2. Log in through CILogon with your NCSA username, password, and Duo MFA.

3. Open the Interactive Apps menu and click Desktop.

4. Configure the job settings and click Launch:
   - Container image: keep default
   - Account: `beyi-delta-gpu`
   - Partition: `GPUA40x4-interactive`
   - Duration: `00-01:00:00`
   - Reservation: leave empty if none
   - CPUs: `2`
   - RAM: `16GB`
   - GPUs: `1`

5. Wait for the job status to change from "starting" to "running" in My Interactive Sessions. 
![starting](https://docs.ncsa.illinois.edu/systems/delta/en/latest/_images/desktop-starting.png)
![running](https://docs.ncsa.illinois.edu/systems/delta/en/latest/_images/desktop-connect.png)
## 2. Load VMD Module
Click "Connect to Desktop" to access the Linux graphical interface. Open a terminal and run:

```bash

module load 
source activate /projects/beyi/sw/conda/envs/vmdplugin
module load vmd/1.9.4lm
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
vmd
```

## 3. Load Trajectory and Visualization State

1. Load the `MinCell_1.lm` trajectory file
2. Load the visualization state `render.vmd`

Explore the trajectory visualization as needed. 

![cell_traj_snapshot](./figures/VMD_render.png)