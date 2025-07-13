# Simulating DNA using btree_chromo and LAMMPS

## Description:
<img align="right" width="300" src="./figures/1. Introduction to simulation with btree_chromo and LAMMPS/initial_state.png">

We will walk you through how to set up and run a simulation using the program btree_chromo, on the Delta HPC cluster. We will simulate the DNA dynamics of the Minimal Cell, including DNA replication, disentanglement of daughter chromosomes, and partitioning of daughter chromosomes into their respective daughter volumes. The coarse-grained model of the DNA, ribosomes and cell membrane will be discussed, as well as the use of LAMMPS to perform energy minimizations and Brownian dynamics. We will also go into greater detail about how we model biological mechanisms such as SMC looping and topoisomerase. You will get a chance to visualize and analyze a simulation trajectory in VMD.

*This tutorial was prepared for the 2nd edition of the STC QCB Summer School, held during July of 2025.*

## Outline of tutorial:

1. Introduction to DNA simulation with btree_chromo and LAMMPS
2. Setting up and running your first simulation on Delta
3. Modeling DNA initial structure and replication
4. Modeling chromosome dynamics
5. Understanding btree_chromo commands
6. Visualization with VMD

Most of the content of this tutorial, especially the implmentation of energy terms for the DNA polymer, DNA disentanglement, and general procedure for simulating Brownian dynamics and energy minimization with LAMMPS on a GPU, has been submitted as part of a manuscript[^thornburg2025] which is currently under review. The content on SMC blocking/bypassing and daughter chromosome partitioning without the need for an additional fictitious force is a work in progress. 

## 1. Introduction to DNA Simulation with btree_chromo and LAMMPS

Here, we simulate DNA replication and dynamics using the C++ program btree_chromo, available online at  [github.com/brg4/btree_chromo](https://github.com/brg4/btree_chromo). This program was created mainly for the purposes of simulating the minimal cell chromosome, but it can be used to simulate any circular chromosome. The main purpose of the program is to model replication states of the chromosome, as well as perform simulation of chromosome dynamics by calling [LAMMPS](https://www.lammps.org/#gsc.tab=0) (Large-scale Atomic/Molecular Massively Parallel Simulator), a molecular dynamics program from Sandia National Laboratories.  

<img align="right" width="300" src="./figures/1. Introduction to simulation with btree_chromo and LAMMPS/DNA_model_0.png">

The DNA that btree_chromo simulates is coarse-grained at a 10 bp resolution. This means that a single, 3.4 nm diameter bead is used to represent 10 base pairs. We also use beads to represent ribosomes (10 nm), and the cell membrane. The program emulates the effects of SMC (structural maintenence of chromosomes) proteins that extrude loops of DNA to effect chromosome organization, as well as type II topoisomerases which allow for DNA strand crossings when they become tangled.

ZLS group member Ben Gilbert (recently graduated) wrote the program and used it to investigate chromosome disentanglement of replicated DNA strands and partitioning of the daughter strands into separate halves of a toy model cell (50 kbp genome), as well as create chromosome contact maps based on the simulated trajectories[^gilbert2023]. Ongoing research using btree_chromo involves integrating btree_chromo with Lattice Microbes, as well as investigating disentanglement and partitioning of a full minimal cell model (543 kbp genome). 

<div align="center">

| <img src="./figures/1.%20Introduction%20to%20simulation%20with%20btree_chromo%20and%20LAMMPS/division_2chromo_5000bp_separation.png" width="300"/>  | <img src="./figures/1.%20Introduction%20to%20simulation%20with%20btree_chromo%20and%20LAMMPS/partition.png" width="300"/> |
|:--:|:--:|
| Figure 1: Cell division of 50 kbp toy model. A repulsive force between the two circular DNA strands (blue, red) helps to partition DNA into the two sides of the cell, accompanying the cell membrane (green) shape change during cell division. | Figure 2: DNA partitioning of 543 kbp model. SMC looping and topoisomerase action has been performed corresponding to the biological time of DNA replication (~60 mins), along with Brownian dynamics timesteps that amount to 20 ms of biological time |

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

This bash script creates a workspace `/projects/bddt/${USER}/btree_chromo_workspace`. It also copies the `examples` directory into the workspace, which contains example input and output files for **btree_chromo**.

**Running btree_chromo via non-interactive job**

Alternatively, we can run our btree_chromo simulation via a bash script that submits a job to Delta such that it runs non-interactively. This has the advantage over the interactive session in that it will perform all the same commands above, while not putting you into the Apptainer shell (the command-line interface within the container). 

You will still obtain output in the folder `/projects/bddt/${USER}/btree_chromo_workspace` in the file `full_model.lammpstrj` but you will no longer get the continous output stream from btree_chromo clogging up your terminal.

Run the command
```bash
sbatch /projects/bddt/${USER}/btree_chromo_workspace/run_btree_chromo.sh
```
If you get `sbatch: error: Unable to open file /projects/bddt/${USER}/btree_chromo_workspace/run_btree_chromo.sh`, please do

```bash
cp /projects/bddt/DNA/files/run_btree_chromo.sh /projects/bddt/${USER}/btree_chromo_workspace/run_btree_chromo.sh
```
and then you should be able to run the `sbatch` command.

If you would like to monitor the progress of your job, you can do the command

```bash
squeue -u ${USER}
```
which will shoud you your jobid, how long the job has been running, as well as what partition and GPU node it is running on.

You can also do 
```bash
ls -lh /projects/bddt/${USER}/btree_chromo_workspace/
```
to see when the last time your `full_model.lammpstrj` was updated, as well as its size. 


> [!TIP]
RNG seed...
>

> [!TIP]
Job run time...
> 


## 3. Modeling DNA initial structure and replication 
In this section you will learn how we generate initial conditions for the chromsome and ribosome bead positions. You will also learn about DNA replication in the minimal cell, and how to represent replication states with a binary tree model.

### Growing the DNA
<img align="right" width="300" src="./figures/3. Modeling the minimal cell/sc_growth_composite_0.png">

At the start of every simulation, we need initial configuration, i.e. coordinates for the DNA, ribosomes, and cell membrane. Initial configurations for the chromosome are generated using a midpoint-displacement algorithm that creates three-dimensional, closed curves formed from overlapping spherocylinder segments. We assume a spherical cell with a known ribosome distribution (nearly randomly distributed, according to Cryo-ET), and "grow in" a self-avoiding chain of these spherocylinders. This process involves adding segments iteratively while avoiding overlaps with ribosomes and preventing knots. Spherical monomers are then interpolated along the spherocylinders. This method accurately models the "fractal globule" chromosome configuration present in Syn3A cells.

See the figure on the right for  a schematic of algorithm used to generate initial conditions of the chromosome. The beads coordinates generated by this algorithm are somewhat jagged, but an energy minimization will relax the structure. We will discuss energy minimizations in the next section. You can double click the image to open up a larger version in a different tab.

The first part of the python script we ran above generated initial configurations for the DNA and ribosomes. The code for performing this algorithm is available at [github.com/brg4/sc_chain_generation](https://github.com/brg4/sc_chain_generation). For the DNA in our simulations we will run below, we used coordinates generated from this program. If one would like to generate coordinates for DNA, as well as ribosomes, one should download, compile and run `sc_chain_generation` from the github link above. In the following simulations, ribosomes have not been included, but it is straightforward to include those: just vary the number of obstacles `N_o` in the input script (`.inp` file) for `sc_chain_generation`. For your convenience, an input script for `sc_chain_generation` for generating a 54338 bead chromosome with 500 ribosomes has been included in this github repository, in `SummerSchool_2024/DNA/files/Syn3A_chromosome_init.inp`. The program `sc_chain_generation` can output the coordinates in either a `.bin`, `.dat`, or `.xyz` file format, the first of which is meant to be read by btree_chromo, and the last of which is human readable and easily read by VMD. To run the input file, you can do  `/path/to/sc_chain_generation/src/gen_sc_chain --i_f=${input_fname} --o_d=${outputDirectory} --o_l=Syn3A_chromosome_init --s=10 --l=${log_fname} --n_t=8 --bin --xyz`.

### Modeling Replication States

The JCVI-syn3A minimal cell has a 543379 bp (543 kbp) genome comprised of 493 genes. This means an unreplicated chromosome is represented as a circular polymer of 54338 beads. Replication begins at a location on the genome called the origin (_Ori_), proceeds along the DNA in the clockwise and counterclockwise directions with Y-shaped structures (Fork), and ends at the terminal site, also called the terminus (_Ter_).  It turns out the replication states of the minimal cell aren't that interesting: it undergoes one replication initiation event per cell cycle, which means it starts with one unreplicated circular chromosome, and replication proceeds from _Ori_ to _Ter_ until we have two complete circular chromosomes.

<img align="center" width="600" src="./figures/3. Modeling the minimal cell/topology_simple.png">

**Figure 3: Representing replication states.**  _Ori_, _Ter_, and Forks given in red, orange, and magenta respectively.

We can use a binary tree to represent the replication state of the chromosome. An unreplicated chromosome is represented by a single node, which we call the mother (m). When unreplicated, the chromosome is circular. As replication proceeds (starting from the _Ori_, along the Forks), the mother branches into two nodes, which we call the left and right daughters (ml and mr). The structure is now no longer circular; it is now called a "theta structure" due to its resemblance to the Greek letter $\theta$. Both the left and right daughters have their own _Ori_'s, so in principle, they could begin to replicate too. 

In the figure below we illustrate replication of a 100 bead chromsome. The origin is depicted in red, the terminus, in orange, and forks in magenta. For each of the four stages of replication, we show the theta structure, the binary tree representation, the physical structure, and finally the bond topology. The bond topology displays all monomers using a colorbar, with the orange arc representing the bond at the terminus, and the pink arcs representing the bonds between the strand of newly created beads at the forks.

<img align="center" width="1000" src="./figures/3. Modeling the minimal cell/replication_topology_0.png">

**Figure 4: Different ways of representing replication states.**  For an unreplicated (left) and partially replicated (right) chromosome, the theta structure is shown in the top-left, the binary tree representation in the top-middle, the physical model in the top-right, and the bond topology of the physical model in the bottom. _Ori_, _Ter_, and Forks given in red, orange, and magenta respectively.

For our simulations, we implement the "train-track" model of bacterial DNA replication[^gogou2021], where replisomes independently move along the opposite arms of the mother chromosome at each replication fork, replicating the DNA. There is another model called the "replication factory" model, but since Syn3A has so few regulatory mechanisms, this second one unlikely. (Plus, the train track model is also more consistent with our understanding of replication initiation[^thornburg2022].) In our implementation, new monomers are added to the left and right daughter chromosomes during replication by creating pairs of monomers centered around the corresponding position of the mother chromosome's monomers. 

<img align="center" width="800" src="./figures/3. Modeling the minimal cell/replication_topology_partB_horizontal_rep_only_0.png">

**Figure 5: Replication with train-track model.**  Starting with an unreplicated Syn3A chromosome (543,379 bp) inside a 200 nm radius cell with 500 ribosomes (not shown), 20,000 bp (2,000 monomers) were replicated using the train-track model (refer to the schematic). Circles are used to highlight the origins of replication (Oris), termination sites (Ter), and replication forks in the replicated system.

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

Both Langevin and Brownian dynamics can be used to correctly sample the NVT ensemble, but Brownian dynamics is preferred in our case since it allows us to take comparatively large time steps. Brownian dynamics is also sometimes called overdamped Langevin dynamics. This approximation is valid for timesteps that satisfy $\Delta t \gg m_i/\gamma_i$.

<img align="right" width=250 src="./figures/4. Modeling chromosome dynamics/bdry_alt.svg">

### SMC looping and topoisomerases

During the genome reduction process of Syn3A, guided by transposon mutagenesis studies on the original JCVI-syn1.0 genome and its intermediate reduced versions, it was found that structural maintenence of chromosomes (SMC) proteins were essential. However, the effect of SMC looping during the minimal cell replication cycle is not fully understood. While magnetic tweezer experiments have been done to determine loop extrusion step size of ~200 bp[^ryu2022], and simulations indicate an extrusion frequency of ~2.5 steps/s[^nomidis2022], we have limited experimental results for SMC looping in the crowded in-cell environment. 

<img align="center" height=300 src="./figures/4. Modeling chromosome dynamics/ganji.png">

**Figure 6: Real-time imaging of DNA loop extrusion by SMC complex.**  A series of snapshots shows DNA loop extrusion intermediates cuased by an SMC dimer on a SxO-stained double-tethered DNA strand. A constant flow at a large angle to DNA axis stretches extruded loop and maintains DNA in imaging plane. Adapted from Ganji et al[^ganji2018].

<img align="right" width=250 src="./figures/4. Modeling chromosome dynamics/looping_bidirectional.png">

Ganji et al. directly visualized the process by which condensin (aka, an SMC dimer) complexes extrude DNA into loops[^ganji2018]. They demonstrated that a single condensin can pull in DNA from one side at a force-dependent rate, supporting the loop extrusion model as a mechanism for chromosome organization. This finding provides strong evidence that SMC protein complexes like condensin actively shape the spatial arrangement of the genome.

The simulation methodology we use for SMC looping is that of Bonato and Michieletto, in which DNA loops are created by adding harmonic bonds bewteen "anchor" and "hinge" monomers[^bonato2021]. Updating the hinge locations causes loop extrusion, as illustrated in the figure to the right. In the simulation, a bead is considered a candidate for the hinge only if it is within a 50 nm radius of the anchor bead; this is because the SMC dimer has a long axis of ~50 nm, and thus can only physically couple DNA beads a maximum 50 nm apart[^nomidis2022].

| Parameter | Description |
| --- | --- |
| Total number of loops | Number of active anchor+hinge pairs that are extruding loops. We know that there are ~100 SMC dimers[^gilbert2023]. |
| Loop extrusion frequency (s^-1) | How often does loop extrusion occur? Our best estimate is around every 0.4 s[^nomidis2022]. |
| Unbind/Rebind frequency (s^-1) | How often does the anchor move to a new location? |
| Extrusion step size (bp) | ~200 bp[^ryu2022]|


Also found to be essential were topoisomerases. There is evidence for coordination between topoisomerases and SMC complexes[^zawadzki2015]. For our simulations, topoisomerase is modeled by periodically running a set of minimizations and Brownian dynamics steps with DNA-DNA pair interactions replaced by soft potentials, which permits strand-crossings.

<img align="right" width=250 src="./figures/4. Modeling chromosome dynamics/topo.png">

<img align="right" width=250 src="./figures/4. Modeling chromosome dynamics/topo2.png">

<img align="right" width=250 src="./figures/4. Modeling chromosome dynamics/SMC_block_bypass.png">

## 5. Understanding btree_chromo Commands

Need to update this section

## 6. Visualization and analysis with VMD

You will now copy over the .lammpstrj files from Delta to your local machine in order to visualize them in vmd. Open up a new terminal, which we will call **Local terminal**.

We can use VMD on DELTA

In the **Local terminal**:

```bash
scp $USERNAME@login.delta.ncsa.illinois.edu:/projects/bddt/$USERNAME/btree_chromo_workspace/\\*.lammpstrj .

```

(Alternatively, instead of `.` you may choose to specify path to a local directory.)
We will also copy over the .tcl scripts which will create nice representations for the DNA and boundary particles.

In the **Local terminal**:

```bash
scp $USERNAME@login.delta.ncsa.illinois.edu:/projects/bddt/$USERNAME/btree_chromo_workspace/\\*.tcl .

```

In VMD, open the VMD TkConsole and do (Extensions->Tk Console). In the Tk Console, do `source full_model.tcl`. 

**Important considerations for Windows Users:**
For those using a Windows machine, you will need to make sure your environmental variables for LAMMPS are set correctly by setting them via command line outside of VMD. Before starting VMD, in Windows Command Shell, please do the following:
```bash
setx LAMMPSDUMMYPOS "$xd,$yd,$zd"
setx LAMMPSMAXATOMS "200000"
setx LAMMPSREMAPFIELDS "vx=c_id_track,vy=c_type_track"
```
Each entry should produce “SUCCESS: Specified value was saved.”

This workaround is only needed on Windows VMD (i.e. not on Linux and Mac VMD). This issue will be addressed in upcoming VMD releases. As of writing this, the latest VMD is Version 1.9.4.

VMD Render stuff

| Monomer type | Color | Bead Size |
| --- | --- | --- |
| DNA | viridis color scale | 13.0 |
| Ori | red | 39.0 |
| Ter | orange | 39.0 |
| Fork | magenta | 39.0 |
| Ribosome | mauve | 70.0 |
| Boundary | gray | 32.5 |
| Anchor | black | 19.5 |
| Hinge | white | 19.5 |

## References
[^gilbert2023]: Gilbert, Benjamin R., Zane R. Thornburg, Troy A. Brier, Jan A. Stevens, Fabian Grünewald, John E. Stone, Siewert J. Marrink, and Zaida Luthey-Schulten. “Dynamics of Chromosome Organization in a Minimal Bacterial Cell.” Frontiers in Cell and Developmental Biology 11 (August 9, 2023). https://doi.org/10.3389/fcell.2023.1214962.
[^gogou2021]: Gogou, Christos, Aleksandre Japaridze, and Cees Dekker. “Mechanisms for Chromosome Segregation in Bacteria.” Frontiers in Microbiology 12 (June 2021). https://doi.org/10.3389/fmicb.2021.685687.
[^thornburg2022]: Thornburg, Zane R., David M. Bianchi, Troy A. Brier, Benjamin R. Gilbert, Tyler M. Earnest, Marcelo C. R. Melo, Nataliya Safronova, et al. “Fundamental Behaviors Emerge from Simulations of a Living Minimal Cell.” Cell 185, no. 2 (January 20, 2022): 345-360.e28. https://doi.org/10.1016/j.cell.2021.12.025.
[^thornburg2025]: Thornburg, Zane R., Andrew Maytin, Jiwoong Kwon, Troy A. Brier, Benjamin R. Gilbert, Enguang Fu, Yang-Le Gao, Jordan Quenneville, Tianyu Wu, Henry Li, Talia Long, Weria Pezeshkian, Lijie Sun, John I. Glass, Angad Mehta, Taekjip Ha, and Zaida Luthey-Schulten. “Bringing the Genetically Minimal Cell to Life on a Computer in 4D.” bioRxiv, June 10, 2025. https://doi.org/10.1101/2025.06.10.658899.
[^ryu2022]: Ryu, Je-Kyung, Sang-Hyun Rah, Richard Janissen, Jacob W J Kerssemakers, Andrea Bonato, Davide Michieletto, and Cees Dekker. “Condensin Extrudes DNA Loops in Steps up to Hundreds of Base Pairs That Are Generated by ATP Binding Events.” Nucleic Acids Research 50, no. 2 (January 25, 2022): 820–32. https://doi.org/10.1093/nar/gkab1268.
[^nomidis2022]: Nomidis, Stefanos K, Enrico Carlon, Stephan Gruber, and John F Marko. “DNA Tension-Modulated Translocation and Loop Extrusion by SMC Complexes Revealed by Molecular Dynamics Simulations.” Nucleic Acids Research 50, no. 9 (May 20, 2022): 4974–87. https://doi.org/10.1093/nar/gkac268.
[^ganji2018]: Ganji, Mahipal, Indra A. Shaltiel, Shveta Bisht, Eugene Kim, Ana Kalichava, Christian H. Haering, and Cees Dekker. “Real-Time Imaging of DNA Loop Extrusion by Condensin.” Science 360, no. 6384 (April 2018): 102–5. https://doi.org/10.1126/science.aar7831.
[^bonato2021]: Bonato, Andrea, and Davide Michieletto. “Three-Dimensional Loop Extrusion.” Biophysical Journal 120, no. 24 (December 2021): 5544–52. https://doi.org/10.1016/j.bpj.2021.11.015.
[^zawadzki2015]: Zawadzki, Pawel, Mathew Stracy, Katarzyna Ginda, Katarzyna Zawadzka, Christian Lesterlin, Achillefs N. Kapanidis, and David J. Sherratt. “The Localization and Action of Topoisomerase IV in Escherichia Coli Chromosome Segregation Is Coordinated by the SMC Complex, MukBEF.” Cell Reports 13, no. 11 (December 22, 2015): 2587–96. https://doi.org/10.1016/j.celrep.2015.11.034.
