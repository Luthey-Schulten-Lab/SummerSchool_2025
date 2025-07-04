# CME-ODE Whole-Cell Model of a Genetically Minimal Cell, JCVI-Syn3A

## Outline:

1. Run CME/ODE Whole-Cell Model in Parallel
2. Minimal Genome and Genetic Information Processes
3. Essential Metabolism
4. Hybrid CME-ODE Algorithm
5. Analysis
6. Discussion
  
## 1. Run CME/ODE Whole-Cell Model in Parallel

We launch independent cell replicates to sample the statistically significant cellular dynamics. To increase the simulation speed, running independent replicates in parallel is inevitable. In our current implementation, we use `mpirun` module to launch simulation python scripts in parallel.  The
main script is `WCM CMEODE Hook.py` calling multiple other python scripts to construct and simulate the genetic information processes, metabolism and their interactions over the whole cell cycle. Essentially, we use LM and build-in Gillepsie algorithm to construct and then simulate the genetic information processes. The ODEs are written/constructed by odecell, a package developed by former group member and simulated by well-known LSODA algorithm in Scipy library.

<p align="center">
  <img src="../figs/figs_WCM/parallel.png" width="900" alt="mpirun to launch parallel simulation"> <br>
  <b>Figure 1. Launching Independent Whole Cell Simulations in Parallel using mpirun module. The user input will be passed to claim the replicates number, time length, and hook interval.</b> <br>
  <b>Each simulation is independent with each other. In current simulation we have four main input files and the trajectories will be stored in CSV files with an additional log file.</b>
</p>

The spatially homogeneous simulations can be efficiently parallelized across up to 25 indepedent cell replicates or more, with each replicate requiring less than 2GB of RAM. On systems equipped with AMD EPYC 7763 “Milan” processors on **[Delta](https://docs.ncsa.illinois.edu/systems/delta/en/latest/index.html)** or Intel Xeon Gold 6154 CPUs @ 3.00 GHz on normal workstation, the parallel simulations of 2 biological hours with communication step of 1 s typically complete within **6 physical hours**.

Due to the time limitation, you will launch 2 minutes simulation of 4 cell replicates. You are encouraged to change the parameters to run longer of more replicates.

### Launch Four Cell Replicates

+ **First**: Open **another** new terminal and login to Delta again

+ **Second**: Navigate and Submit the bash file

    ``` bash
    cd /projects/bddt/$USER/LM/CME/WholeCellModel/programs
    ```

    ```bash
    sbatch mpirun.sh
    ```
  
+ **Third**: Check the status of your job.  
    ```bash
    squeue -u $USER
    ```
    *PD* means waiting to run, *R* running.

    Go to output folder `output_4replicates` 
    
    ``` bash
    cd /projects/bddt/$USER/LM/CME/WholeCellModel/output_4replicates
    ```

Each simulation replicate with index *i* generates:

- `counts_i.csv`: Species count trajectories of metabolites from ODE and genetic particles from CME (units: molecules).
- `SA_i.csv`: Surface area (nm$^2$ or m$^2$) and volume (L) trajectories.
- `Flux_i.csv`: Fluxes through ODE reactions (units: mM/s).
- `log_i.txt`: Log file with timestamps, printed reactions, run times, and any warnings/errors.

Output files are saved to directories defined and created in `mpirun.sh`. Typical CSV file size ranges from **100–200 MB** for a 7200 s simulation with 1 s hook intervals.
### Scripts

<details>
<summary><strong>Click to expand: Derivation of the RDME</strong></summary>

Here is a detailed derivation of the Reaction-Diffusion Master Equation (RDME).

Start from the chemical master equation (CME):
$$
\frac{dP(\mathbf{x},t)}{dt} = \cdots
$$

Then extend it by adding diffusion terms:
$$
\cdots
$$

</details>

### Input Files

Four main input files are used in whole-cell simulation: one Genbank file, one SBML file, and two Excel files. The `syn3A.gb` Genebank file contains the sequences and functions of genes, RNAs, and proteins, and the `Syn3A_updated.xml` SBML file contains the metabolic reactions (reactants, stoichiometries). `The intitial_concentration.xlsx` file contains the initial count/concentrations of proteins and metabolite while `kinetic params.xlsx` contains the kinetic parameters of the GIP and metabolic reactions.

- `syn3A.gb` — GenBank file of JCVI-Syn3A 
  - Obtained from ([NCBI Accession: CP016816](https://www.ncbi.nlm.nih.gov/nuccore/CP016816)).
  - Encodes genome sequence, segmentation, and gene annotations.

- `Syn3A_updated.xml` — Includes metabolites, compartments, reactions, and gene associations.
  - Updated SBML model from [*eLife*, 2019](https://elifesciences.org/articles/36842).
  - Standard file format in system biology.
  - Read in for Reaction and stoichiometries by `rxns_ODE.py`.

- `initial_concentration.xlsx` — Provides initial conditions for proteins, medium, and metabolites.
  -  Update from [*Cell*, 2022](https://www.sciencedirect.com/science/article/pii/S0092867421014884?via%3Dihub#da0010).
  - **Sheet breakdown**:
  - **Comparative Proteomics**: Protein initial counts for CME.
  - **Simulation Medium**: Medium composition for simulation.
  - **Intracellular Metabolites**: Metabolite concentrations in cytoplasm for ODE.
  - **mRNA count**: Initial average count of mRNAs for CME.
  - **Protein Metabolites**: Protein metabolite IDs used in `initiation.py` and `rxns_ODE.py`.

- `kinetic_params.xlsx` — Contains kinetic parameters for ODE reactions, tRNA charging, and gene expression
  - Update from [*Cell*, 2022](https://www.sciencedirect.com/science/article/pii/S0092867421014884?via%3Dihub#da0010).
  - **Sheet Breakdown:**
  - **Central, Nucleotide, Lipid, Cofactor, Transport**: Random binding + convenience rate law reactions (`rxns_ODE.py`).
  - **Non-Random-binding Reactions**: Passive transport and serial phosphorelay reactions (`rxns_ODE.py`).
  - **tRNA Charging**: Aminoacylation parameters (`rxns_CME.py`).
  - **Gene Expression**: Parameters in gene expression (`rxns_CME.py`, `replication.py`, and `GIP_rates.py`)
  - **SSU Assembly, LSU Assembly**: Reactions and rates of SSU and LSU assembly

- `complex_formation.xlsx` — Defines the composition and initial counts of protein complexes.




## 2. Minimal Genome and Genetic Information Processes

## 3. Essential Metabolism

## 4. Hybrid CME-ODE Algorithm

## 5. Analysis

### Run pkl.ipynb then plotting.ipynb on Jupyter Notebook Webpage
+ **First**: Navigate to `./LM/CME/WholeCellModel/analysis` in Jupyter Notebook Webpage.
  
+ **Second**: Open *pkl.ipynb* and **Run All** to serialize the CSV files into one single python pickle file.  
  10 CSV files are generated beforehand for the analysis.
  
+ **Thrid**: Open *plotting.ipynb* and **Run All** to generate the plots.  
  Compare the plots with figures in *Summer_School_CME_2024.pdf*.

## 6. Discussion
