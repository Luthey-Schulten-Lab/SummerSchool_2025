# Simulating DNA using btree_chromo and LAMMPS

## Description:
<img align="right" width="300" src="./figures/1. Introduction to simulation with btree_chromo and LAMMPS/initial_state.png">

The purpose of this tutorial is to give you a crash course on modeling and simulating the chromosome of the minimal cell. 

We will walk you through how to set up and run a simulation using the program btree_chromo, on the Delta HPC cluster. The coarse-grained model of the DNA, ribosomes and cell membrane will be discussed, as well as the use of LAMMPS to perform energy minimizations and Brownian dynamics. We will also go into greater detail about how we model biological mechanisms such as SMC looping and topoisomerase. Finally, you will get a chance to visualize and analyze a simulation trajectory in VMD.

*This tutorial was prepared for the STC QCB Summer School, held during August of 2024.*

## Outline of tutorial:

1. Introduction to DNA simulation with btree_chromo and LAMMPS
2. Setting up and running your first simulation on Delta
3. Modeling the minimal cell chromsome replication and initial structure
4. Modeling chromosome dynamics
5. Understanding btree_chromo commands
6. Visualization and analysis with VMD

## 1. Introduction to DNA Simulation with btree_chromo and LAMMPS

Here, we simulate DNA replication and dynamics using the C++ program btree_chromo, available online at  [github.com/brg4/btree_chromo](https://github.com/brg4/btree_chromo). This program was created mainly for the purposes of simulating the minimal cell chromosome, but it can be used to simulate any circular chromosome. The main purpose of the program is to model replication states of the chromosome, as well as perform simulation of chromosome dynamics by calling [LAMMPS](https://www.lammps.org/#gsc.tab=0) (Large-scale Atomic/Molecular Massively Parallel Simulator), a molecular dynamics program from Sandia National Laboratories.  

<img align="right" width="300" src="./figures/1. Introduction to simulation with btree_chromo and LAMMPS/DNA_model_0.png">

The DNA that btree_chromo simulates is coarse-grained at a 10 bp resolution. This means that a single, 3.4 nm diameter bead is used to represent 10 base pairs. We also use beads to represent ribosomes (10 nm), and the cell membrane. The program emulates the effects of SMC (structural maintenence of chromosomes) proteins that extrude loops of DNA to effect chromosome organization, as well as type II topoisomerases which allow for DNA strand crossings when they become tangled.

ZLS group member Ben Gilbert (recently graduated) wrote the program and used it to investigate chromosome disentanglement of replicated DNA strands and partitioning of the daughter strands into separate halves of a toy model cell (50 kbp genome), as well as create chromosome contact maps based on the simulated trajectories[^gilbert2023]. Ongoing research using btree_chromo involves integrating btree_chromo with Lattice Microbes, as well as investigating disentanglement and partitioning of a full minimal cell model (543 kbp genome). 

<div align="center">

| <img src="./figures/1.%20Introduction%20to%20simulation%20with%20btree_chromo%20and%20LAMMPS/division_2chromo_5000bp_separation.png" width="300"/>  | <img src="./figures/1.%20Introduction%20to%20simulation%20with%20btree_chromo%20and%20LAMMPS/partition.png" width="300"/> |
|:--:|:--:|
| Figure 1: Cell division of 50 kbp toy model. A repulsive force between the two circular DNA strands (blue, red) helps to partition DNA into the two sides of the cell, accompanying the cell membrane (green) shape change during cell division. | Figure 2: DNA partitioning of 543 kbp model. SMC looping and topoisomerase action has been performed corresponding to the biological time DNA replication (~60 mins), along with Brownian dynamics timesteps that amount to 20 ms of biological time |

</div>
Today you will run a simulation using a variant of LAMMPS which utilizes the GPUs on the Delta HPC cluster. We will simulate the entire minimal cell including the effects of SMC proteins, topoisomerase, and Brownian dynamics. 

## 2. Setting up and running your first simulation on Delta
In this section, we will log on to Delta and launch a container which has btree_chromo and LAMMPS already installed. Then, we will start running a simulation of the minimal cell chromosome. 

> [!NOTE]
The reason we are doing this first, is so that the simulation will be left running throughout the tutorial. At the end of this tutorial, we will visualize and analyze the results of our simulations.
> 

**Step 1: Log in to Delta**

Connect to the Delta HPC using the following command. You will need to type your password and do 2FA.

```bash
ssh $USERNAME@login.delta.ncsa.illinois.edu
```

> [!WARNING]
You will need to replace $USERNAME with your own.
> 

**Step 2: Create workspace and copy examples folder**

```bash
bash /projects/bddt/DNA/files/prelaunch_btree_chromo.sh

```

This bash script creates a workspace `/projects/bcuj/${USER}/btree_chromo_workspace`. It also copies the `examples` directory into the workspace, which contains example input and output files for **btree_chromo**.

**Step 3: Launch the container**

```bash
bash /projects/bddt/DNA/files/launch_btree_chromo.sh

```

This will run the container in interactive mode. You should now see the `Apptainer>` prompt which indicates your have entered the container. We will be running btree_chromo and viewing the terminal output in the container.

Do the command

```bash
export LD_LIBRARY_PATH="/usr/local/lib64:/usr/local/lib:/usr/local/nvidia/lib:/usr/local/nvidia/lib64:/.singularity.d/libs"
```
This ensures the program uses the nvidia drivers appropriate for the NVIDIA A100 GPUs on Delta.

**Step 4: Run btree_chromo**

In the container, please navigate into the folder where the btree_chromo executable is located:

```bash
cd /Software/btree_chromo/build/apps

```

Next, run 'full_model.inp'. Do the command:

```bash
./btree_chromo /mnt/full_model.inp

```
This command starts a simulation of the full 543 kbp chromosome in a 500 nm diameter spherical cell. By the end of this tutorial, our goal is to understand what is going on in the simulation, the contents of the input file `full_model.inp`, as well as some background biological knowledge.

## 3. Modeling the Minimal Cell replication and initial structure
In this section you will learn about DNA replication in the minimal cell, and how to represent replication states with a binary tree model. You will also learn how we generate initial conditions for the chromsome and ribosome bead positions.

### Modeling Replication States

The JCVI-syn3A minimal cell has a 543379 bp (543 kbp) genome comprised of 493 genes. This means an unreplicated chromosome is represented as a circular polymer of 54338 beads. DNA replication, even in the minimal cell where the replication machinery retains only the essential components, is still rather complicated. For today, all we need to understand is that replication begins at a location on the genome called the origin (_Ori_), proceeds along the DNA in the clockwise and counterclockwise directions with Y-shaped structures (Fork), and ends at the terminal site, also called the terminus (_Ter_). 

<img align="center" width="600" src="./figures/3. Modeling the minimal cell/topology_simple.png">

**Figure 3: Representing replication states.**  _Ori_, _Ter_, and Forks given in red, orange, and magenta respectively.

We can use a binary tree to represent the replication state of the chromosome. An unreplicated chromosome is represented by a single node, which we call the mother (m). When unreplicated, the chromosome is circular. As replication proceeds (starting from the _Ori_, along the Forks), the mother branches into two nodes, which we call the left and right daughters (ml and mr). The structure is now no longer circular; it is now called a "theta structure" due to its resemblance to the Greek letter $\theta$. Both the left and right daughters have their own _Ori_'s, so in principle, they could begin to replicate too. 

In the figure below we illustrate replication of a 100 bead chromsome. The origin is depicted in red, the terminus, in orange, and forks in magenta. For each of the four stages of replication, we show the theta structure, the binary tree representation, the physical structure, and finally the bond topology. The bond topology displays all monomers using a colorbar, with the orange arc representing the bond at the terminus, and the pink arcs representing the bonds between the strand of newly created beads at the forks.

<img align="center" width="1000" src="./figures/3. Modeling the minimal cell/replication_topology_0.png">

**Figure 4: Different ways of representing replication states.**  For an unreplicated (left) and partially replicated (right) chromosome, the theta structure is shown in the top-left, the binary tree representation in the top-middle, the physical model in the top-right, and the bond topology of the physical model in the bottom. _Ori_, _Ter_, and Forks given in red, orange, and magenta respectively.

Quickly, let's see how btree_chromo handles replication states. In a new terminal log on to Delta again, and navigate to `/projects/bddt/$USERNAME/btree_chromo_workspace/examples/preparing_chromosome`:

```bash
cd /projects/bddt/$USERNAME/btree_chromo_workspace/examples/preparing_chromosome
ls

```

You should see the following files:

| File name | Description |
| --- | --- |
| preparing_chromosome_directives.inp | contains directives (commands to be executed by btree_chromo) |
| chromo_state_0.dat | contains the initial replication state |
| chromo_state_1.dat | contains the final replication state |
| transforms.dat | contains a set of replication events |

Using a text exitor such as **vim**, open 'preparing_chromosome_directives.inp'. Scroll through the file to see all of the commands. Here, **btree_chromo** first creates a chromosome with 1000 monomers, and then applies a series of transforms, printing the replication state at each step, and finally outputs the final state.

> [!WARNING]
Make sure that the "terminate'' at the top is commented to execute all of the commands.
> 

Next, we will run the program in the **Apptainer terminal**. In the **Apptainer terminal** run 'preparing_chromosome_directives.inp'.

In the **Apptainer terminal**:

```bash
./btree_chromo /mnt/examples/preparing_chromosome/preparing_chromosome_directives.inp

```

Check out the terminal output. Notice how **btree_chromo** starts by testing the validity of all the commands and parameters and prints them to the terminal output. It then executes each of the commands. For example, you should see this output corresponding to the first replication transformation:

```bash
COMMAND: transform
	param_0: m_cw100_ccw200

COMMAND: print

printing tree with 2 leaves and 1 forks
fork breakdown: 0 completed, 1 active
leaves:
ml mr
active forks:
m
total_size = 1300
| generation = 0
| rho_t = 300/1000, rho_cw = 100/800, rho_ccw = 200/900
| start = 0, mid = 500, end = 999
| start_link = 999, end_link = 0
| left branch
  | generation = 1
  | start = 0, mid = 500, end = 999
  | start_link = 999, end_link = 0
| right branch
  | generation = 1
  | start = 1000, mid = 1200, end = 1299
  | start_link = 299, end_link = 600

```
This is a very thorough way of keeping track of what each bead corresponds to in your simulation. If one needs to simulate complicated theta structures, then keeping track of replication states this way is helpful. For our purposes, we don't need to worry about the binary tree formalism too much.

For our simulations, we implement the "train-track" model of bacterial DNA replication (Gogou), where replisomes independently move along the opposite arms of the mother chromosome at each replication fork, replicating the DNA. There is another model called the "replication factory" model, but since Syn3A has so few regulatory mechanisms, this second one unlikely. (Plus, the train track model is also more consistent with our understanding of replication initiation[^thornburg2022].) In our implementation, new monomers are added to the left and right daughter chromosomes during replication by creating pairs of monomers centered around the corresponding position of the mother chromosome's monomers. 

<img align="center" width="800" src="./figures/3. Modeling the minimal cell/replication_topology_partB_horizontal_rep_only_0.png">

**Figure 4: Replication with train-track model.**  Starting with an unreplicated Syn3A chromosome (543,379 bp) inside a 200 nm radius cell with 500 ribosomes, 20,000 bp (2,000 monomers) were replicated using the train-track model (refer to the schematic). Circles are used to highlight the origins of replication (Oris), termination sites (Ter), and replication forks in the replicated system.


### Growing the DNA
<img align="right" width="500" src="./figures/3. Modeling the minimal cell/sc_growth_composite_0.png">

At the start of every simulation, we need initial configuration, i.e. coordinates for the DNA, ribosomes, and cell membrane. Initial configurations for the chromosome are generated using a midpoint-displacement algorithm that creates three-dimensional, closed curves formed from overlapping spherocylinder segments. We assume a spherical cell with a known ribosome distribution (nearly randomly distributed, according to Cryo-ET), and "grow in" a self-avoiding chain of these spherocylinders. This process involves adding segments iteratively while avoiding overlaps with ribosomes and preventing knots. Spherical monomers are then interpolated along the spherocylinders. This method accurately models the "fractal globule" chromosome configuration present in Syn3A cells.

We won't use it here, but the code for performing this algorithm is available at [github.com/brg4/sc_chain_generation](https://github.com/brg4/sc_chain_generation). For our simulation, we used coordinates generated from this program. 

See the figure on the right for  a schematic of algorithm used to generate initial conditions of the chromosome. The beads coordinates generated by this algorithm are somewhat jagged, but an energy minimization will relax the structure. We will discuss energy minimizations in the next section.

## 4. Modeling Chromosome Dynamics

The total potential energy for the chromosome/ribosome system is

$$U= \sum_{i=1}^{N_{\mathrm{DNA}}}\left[U_i^b+U_i^t+U_i^a+U_i^s\right] +\sum_{i=1}^{N_{\mathrm{DNA}}-1} \sum_{j=i+1}^{N_{\mathrm{DNA}}} U_{i j}^{\mathrm{DNA}-\mathrm{DNA}}+\sum_{i=1}^{N_{\mathrm{DNA}}} \sum_j^{N_{\text {ribo }}} U_{i j}^{\mathrm{DNA}-\text { ribo }} +\sum_{i=1}^{N_{\text {ribo }}-1} \sum_{j=i+1}^{N_{\text {ribo }}} U_{i j}^{\text {ribo-ribo }} +\sum_{i=1}^{N_{\text {bdry }}} \sum_j^{N_{\mathrm{DNA}}} U_{i j}^{\text {bdry-DNA }}+\sum_{i=1}^{N_{\text {bdry }}} \sum_j^{N_{\text {ribo }}} U_{i j}^{\text {bdry-ribo }}.$$

The energies for the bending, twisting, stretching and excluded volume interactions are shown below.

|**Bending:** |**Twisting** | **Stretching** | **Excluded Volume** |
|:--:|:--:|:--:|:--:|
| <img src="./figures/4.%20Modeling%20chromosome%20dynamics/DNA_model_bending_0.png/" width="300"/>  | <img src="./figures/4.%20Modeling%20chromosome%20dynamics/DNA_model_twisting_0.png/" width="300"/> |  <img src="./figures/4.%20Modeling%20chromosome%20dynamics/DNA_model_stretching_0.png/" width="200"/> |  <img src="./figures/4.%20Modeling%20chromosome%20dynamics/DNA_model_LJ_0.png/" width="120"/> |
|$U_i^b=\kappa_b\left[1-\cos \left(\pi-\theta_i\right)\right]$|$U_i^t=\kappa_t\left[1-\cos \left(\alpha_i+\gamma_i\right)\right]$|$U_i^s= -\frac{\kappa_s L_0^2}{2} \log \left[1-\left(l_i / L_0\right)^2\right]$ $+4 \epsilon_s\left[\left(\frac{\sigma_s^s}{l_i}\right)^{12}-\left(\frac{\sigma_s}{l_i}\right)^6\right]$ $\times \Theta\left(2^{\frac{1}{6}} \sigma_s-l_i\right)$  |$U_{i j}^{e . v .}=  4 \epsilon_{e . v}\left[\left(\frac{\sigma_{e . v}}{r_{i j}}\right)^{12}-\left(\frac{\sigma_{e . v}}{r_{i j}}\right)^6\right]$ $\times  \Theta\left(2^{\frac{1}{6}} \sigma_{{e.v. }}-r_{i j}\right)$|

> [!NOTE]
In the current version of btree_chromo, we neglect the twisting potential for neighboring beads (the corresponding forces are not calculated and are not used during our simulations). The reason for this has to do with how minimizations work in LAMMPS.
> 

In btree_chromo, we simulate dynamics of the DNA and ribosomes using a GPU-accelerated version of the Brownian dynamics integrator, which performs time-integration to update the coordinates of each of the beads. Let's quickly go through how the Brownian dynamics integrator works. First, recall how the acceleration of each bead is related to the net force on that bead and it's mass, which is given by Newton's second law, $F_{\text{net}}=ma$. 

For bead $i$, which has mass $m_i$, there are three forces acting on the bead: the "system" force associated with the interactions between beads, $F_{\text{system}} = -\nabla U_i $, as well as the drag force $F_{\text{drag}}$ and random force $F_{\text{rand}}$. So, Newton's second law reads

<img align="center" height=50 src="./figures/4. Modeling chromosome dynamics/Newton_eq.png">

The drag force is proportional to the bead's velocity, $\mathbf{v} = \text{d}\mathbf{x}/\text{d}t$, as well as its translational damping constant $\gamma_i$, given by the Einstein-Stokes equation: 

<img align="center" height=25 src="./figures/4. Modeling chromosome dynamics/Stokes-Einstein_eq.png">

Plugging in $F_{\text{drag}} = -\gamma_i v_i$ and dividing by $\gamma_i$, we obtain the Langevin equation of motion:

<img align="center" height=50 src="./figures/4. Modeling chromosome dynamics/Langevin_eq.png">

In the large friction limit, i.e. where $\gamma_i$ is very large, we can neglect the left hand side of the Langevin equation, and we obtain the equation for Brownian motion:

<img align="center" height=50 src="./figures/4. Modeling chromosome dynamics/Brownian_eq.png">

LAMMPS uses this equation to update positions, according to a simple first order integration scheme known as the Euler-Maruyama method (basically, the Euler method).

Both Langevin and Brownian dynamics can be used to correctly sample the NVT ensemble, but Brownian dynamics is preferred in our case since it allows us to take comparatively large time steps. Brownian dynamics is also sometimes called overdamped Langevin dynamics. This approiximation is valid for timesteps that satisfy $\Delta t \gg m_i/\gamma_i$.

### SMC looping and topoisomerases

<img align="right" width=300 src="./figures/4. Modeling chromosome dynamics/DNA_model_looping_0.png">

During the genome reduction process of Syn3A, guided by transposon mutagenesis studies on the original JCVI-syn1.0 genome and its intermediate reduced versions, it was found that structural maintenence of chromosomes (SMC) proteins were essential. However, the effect of SMC looping during the minimal cell replication cycle is not fully understood. While magnetic tweezer experiments have been done to determine loop extrusion step size of ~200 bp/s[^ryu2022], and simulations indicate an extrusion frequency of ~2.5 steps/s[^nomidis2022], we have limited experimental results for SMC looping in the crowded in-cell environment. 

The simulation methodology we use for SMC looping is that of Bonato and Michieletto, in which DNA loops are created by adding harmonic bonds bewteen "anchor" and "hinge" monomers. Updating the hinge locations causes loop extrusion, as illustrated in the figure to the right. 

| Parameter | Description |
| --- | --- |
| Total number of loops | Number of active anchor+hinge pairs that are extruding loops. We know that there are ~100 SMC dimers[^gilbert2023]. |
| Loop extrusion frequency (s^-1) | How often does loop extrusion occur? Our best estimate is around every 0.4 s[^nomidis2022]. |
| Unbind/Rebind frequency (s^-1) | How often does the anchor move to a new location? |
| Extrusion step size (bp) | ~200 bp[^ryu2022]|


Also found to be essential were topoisomerases. There is evidence for coordination between topoisomerases and SMC complexes[^zawadzki2015]. For our simulations, topoisomerase is modeled by periodically running a set of minimizations and Brownian dynamics steps with DNA-DNA pair interactions replaced by soft potentials, which permits strand-crossings.

## 5. Understanding btree_chromo Commands

Let's take a look in `full_model.inp` to get a feel for how to write input files for btree_chromo. You can either use a text editor to open the file in your terminal window, or simply view it on github at [files/full_model.inp](https://github.com/enguangfu/SummerSchool_2024/blob/main/DNA/files/full_model.inp).

```bash
switch_skip_runs:F

btree_prng_seed:10
replicator_prng_seed:10
new_chromo:54338

load_BD_lengths:/mnt/in_BD_lengths_LAMMPS_test.txt
load_mono_coords:/mnt/x_chain_Syn3A_chromosome_init_rep00001.bin,row
load_bdry_coords:/mnt/2500A_bdry.bin,row
prepare_simulator:/mnt/logfile0.log
simulator_set_prng_seed:42
simulator_set_nProc:8
simulator_set_DNA_model:/Software/btree_chromo/LAMMPS_DNA_model_kk
simulator_set_output_details:/mnt/,full_model
simulator_set_delta_t:1.0E+5
```
These commands create a new chromosome and load in coordinates for the DNA and boundary, creates output files, sets the LAMMPS random number generator seed, and sets the timestep to 0.1 ns

```bash
switch_twisting_angles: F
```
This command disables twisting angles between DNA monomers (default T). Turning off twisting removes the $U_i^t$ and $U_i^a$ (cosine potentials for twisting and aligning).

```bash
simulator_minimize_soft_harmonic:500
```
This performs an energy minimization according to the [conjugate gradient algorithm](https://docs.lammps.org/minimize.html).

```bash
sync_simulator_and_system
set_initial_state
transform:m_cw1360_ccw1360
set_final_state
map_replication
sys_write_sim_read_LAMMPS_data:/mnt/data.lammps_0
```
These commands copy the LAMMPS simulation bead coordinates into the binary tree data structure in btree_chromo, replicates by moving both forks by 1360 beads, maps the binary tree structure back onto the bead coordinates, and then updates the bead coordinates in the LAMMPS simulation.

```bash
simulator_run_loops:F,100,200000,50000,20000,append,nofirst
```
This command runs Brownian dynamics with the hard/FENE potential for 200k timesteps with 100 loops randomly placed, while printing thermodynamic information every 50k steps and dumping the coordintaes every 20k steps.

If we take a look at loop_params.txt, we find that it specifies the minimum number of monomers separating anchor and hinge, the distribution of extrusion steps, hinge unbinding probability and grab radius, loop update frequency and topoisomerase relaxation frequency.

```bash
# system parameters

# minimum distance separating anchor and hinge on strand - [# monomers]
min_dist=5

# distribution family for steps
family=poisson
# average extrusion distance [# monomers]
ext_avg=20.0
# max extrusion distance [# monomers]
ext_max=30

# hinge unbinding probability
p_unbinding=0.0
# grab radius - [A]
r_g=5.0E+2

# simulator parameters

# loop update frequency - [# timesteps]
freq_loop=10000

# topoisomerase relaxation (soft DNA-DNA pairs) frequency - [# timesteps]
freq_topo=50000
# topoisomerase relaxation (soft DNA-DNA pairs) interval - [# timesteps]
dNt_topo=50000

```

## 6. Visualization and analysis with VMD

You will now copy over the .lammpstrj files from Delta to your local machine in order to visualize them in vmd. Open up a new terminal, which we will call **Local terminal**.

In the **Local terminal**:

```bash
scp $USERNAME@login.delta.ncsa.illinois.edu:/projects/bddt/$USERNAME/btree_chromo_workspace/\\*.lammpstrj .

```

(Alternatively, instead of `.` you may choose to specify path to a local directory.)
We will also copy over the .tcl scripts which will create nice representations for the DNA, ribosomes and boundary particles.

In the **Local terminal**:

```bash
scp $USERNAME@login.delta.ncsa.illinois.edu:/projects/bddt/DNA/files/\\*.tcl .

```

In VMD, delete the previous two molecules. Open the VMD TkConsole and do (Extensions->Tk Console). In the Tk Console, do `source load_example_full.tcl`. See the DNA polymer replicate!

| Monomer type | Color | Bead Size |
| --- | --- | --- |
| DNA | viridis color scale | 13.0 |
| Ori | red | 26.0 |
| Ter | orange | 26.0 |
| Fork | magenta | 26.0 |
| Ribosome | mauve | 70.0 |
| Boundary | gray | 32.5 |
| Anchor | black | 19.5 |
| Hinge | white | 19.5 |

### Calculate and plot the Radius of Gyration
We'll write a small script in the Tcl Console to calculate the radius of gyration for each frame of the trajectory and store the results.

Define Variables for the output file and number of frames:

```bash
set outfile [open "rgyr_vs_frame.txt" w]
set num_frames [molinfo top get numframes]
```

Loop Over All Frames:

```bash
for {set i 0} {$i < $num_frames} {incr i} {
    animate goto $i
    set DNA [atomselect top {vy>2}]
    set rg [measure rgyr $DNA]
    puts $outfile "$i $rg"
}
```

This script calculates the radius of gyration for all atoms in each frame and stores the values in the list rg_values.

Next, save the data to the designated output file:

```bash
close $outfile
```
To visualize the change in RoG over time, we can use the following python script:
```python
import matplotlib.pyplot as plt

# Read the data from the file
frames = []
rgy_values = []

with open("rgyr_vs_frame.txt", "r") as file:
    for line in file:
        frame, rg = line.split()
        frames.append(int(frame))
        rgy_values.append(float(rg))

# Get the initial radius of gyration value for percentage calculation
initial_rg = rgy_values[0]

# Calculate the percentage values
percentage_values = [(rg / initial_rg) * 100 for rg in rgy_values]

# Create the plot
fig, ax1 = plt.subplots(figsize=(10, 6))

# Plot radius of gyration on the primary y-axis
ax1.plot(frames, rgy_values, marker='o', linestyle='-', color='b', label='Radius of Gyration', linewidth=1)
ax1.set_xlabel('Frame')
ax1.set_ylabel('Radius of Gyration (Å)', color='b')
ax1.tick_params(axis='y', labelcolor='b')
ax1.grid(True)

# Create a secondary y-axis for percentage values
ax2 = ax1.twinx()
ax2.plot(frames, percentage_values, marker='o', linestyle='--', color='r', label='Percentage of Initial Value', linewidth=1)
ax2.set_ylabel('Percentage of Initial Radius of Gyration (%)', color='r')
ax2.tick_params(axis='y', labelcolor='r')

# Add titles and labels
plt.title('Radius of Gyration vs. Frame')

# Save the plot as a PNG file
plt.savefig("radius_of_gyration_vs_frame.png")

# Show the plot
plt.show()
```

## References
[^gilbert2023]: Gilbert, Benjamin R., Zane R. Thornburg, Troy A. Brier, Jan A. Stevens, Fabian Grünewald, John E. Stone, Siewert J. Marrink, and Zaida Luthey-Schulten. “Dynamics of Chromosome Organization in a Minimal Bacterial Cell.” Frontiers in Cell and Developmental Biology 11 (August 9, 2023). https://doi.org/10.3389/fcell.2023.1214962.
[^thornburg2022]: Thornburg, Zane R., David M. Bianchi, Troy A. Brier, Benjamin R. Gilbert, Tyler M. Earnest, Marcelo C. R. Melo, Nataliya Safronova, et al. “Fundamental Behaviors Emerge from Simulations of a Living Minimal Cell.” Cell 185, no. 2 (January 20, 2022): 345-360.e28. https://doi.org/10.1016/j.cell.2021.12.025.
[^ryu2022]: Ryu, Je-Kyung, Sang-Hyun Rah, Richard Janissen, Jacob W J Kerssemakers, Andrea Bonato, Davide Michieletto, and Cees Dekker. “Condensin Extrudes DNA Loops in Steps up to Hundreds of Base Pairs That Are Generated by ATP Binding Events.” Nucleic Acids Research 50, no. 2 (January 25, 2022): 820–32. https://doi.org/10.1093/nar/gkab1268.
[^nomidis2022]: Nomidis, Stefanos K, Enrico Carlon, Stephan Gruber, and John F Marko. “DNA Tension-Modulated Translocation and Loop Extrusion by SMC Complexes Revealed by Molecular Dynamics Simulations.” Nucleic Acids Research 50, no. 9 (May 20, 2022): 4974–87. https://doi.org/10.1093/nar/gkac268.
[^bonato2021]: Bonato, Andrea, and Davide Michieletto. “Three-Dimensional Loop Extrusion.” Biophysical Journal 120, no. 24 (December 2021): 5544–52. https://doi.org/10.1016/j.bpj.2021.11.015.
[^zawadzki2015]: Zawadzki, Pawel, Mathew Stracy, Katarzyna Ginda, Katarzyna Zawadzka, Christian Lesterlin, Achillefs N. Kapanidis, and David J. Sherratt. “The Localization and Action of Topoisomerase IV in Escherichia Coli Chromosome Segregation Is Coordinated by the SMC Complex, MukBEF.” Cell Reports 13, no. 11 (December 22, 2015): 2587–96. https://doi.org/10.1016/j.celrep.2015.11.034.
