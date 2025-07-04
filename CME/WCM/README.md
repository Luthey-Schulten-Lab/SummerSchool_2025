# CME-ODE Whole-Cell Model of a Genetically Minimal Cell, JCVI-Syn3A

## Outline:
1. [Run CME/ODE Whole-Cell Model in Parallel](#1-run-cmeode-whole-cell-model-in-parallel)
2. [Minimal Genome and Genetic Information Processes](#2-minimal-genome-and-genetic-information-processes)
3. [Essential Metabolism](#3-essential-metabolism)
4. [Hybrid CME-ODE Algorithm](#4-hybrid-cme-ode-algorithm)
5. [Analysis](#5-analysis)
6. [Discussion](#6-discussion)
  
## 1. Run CME/ODE Whole-Cell Model in Parallel

We launch independent cell replicates to sample the statistically significant cellular dynamics. To increase the simulation speed, running independent replicates in parallel is inevitable. In our current implementation, we use `mpirun` module to launch simulation python scripts in parallel.  The
main script is `WCM CMEODE Hook.py` calling multiple other python scripts to construct and simulate the genetic information processes, metabolism and their interactions over the whole cell cycle. Essentially, we use LM and build-in Gillepsie algorithm to construct and then simulate the genetic information processes. The ODEs are written/constructed by odecell, a package developed by former group member and simulated by well-known LSODA algorithm in Scipy library.

### Scripts

<details>
<summary><strong>Click to EXPAND: Explanation of All Python Scripts </strong></summary>

#### Main Driver

- `WCM_CMEODE_Hook.py` — Main script to launch the CMEODE simulation.

#### Simulation Core

- `species_counts.py` — Python Class for managing species count data in the hook algorithm.
- `integrate.py` — Performs ODE integration using *lsoda* from Scipy.
- `initiation.py` — Initializes constants, time traces, and cell membrane.
- `communicate.py` — Syncs CME and ODE states, computes costs, and updates the membrane.
- `hookSolver_CMEODE.py` — Defines `hookSimulation`, the main hook interval operation.
- `hook_CMEODE.py` — Manages operations that communicate CME and ODE during hooks.
- `filesaving.py` — Exports time traces, surface area, and fluxes to CSV files.

#### Biology Modules

- `rxns_CME.py` — Adds genetic information process (GIP) reactions to CME (e.g., replication, transcription, translation, tRNA charging).
- `cme_complexation.py` — Adds protein complex formation reactions to CME.
- `rxns_ODE.py` — Builds the ODE system using `odecell` (metabolic parameters from `kinetic_params.xlsx`).
- `replication.py` — Defines DNA replication initiation and reactions (GIP parameter from `kinetic_params.xlsx`).
- `GIP_rates.py` — Computes reaction rates for GIP processes (GIP parameter from `kinetic_params.xlsx`).

</details>

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


### Input Files

Four main input files are used in whole-cell simulation: one Genbank file, one SBML file, and two Excel files. The `syn3A.gb` Genebank file contains the sequences and functions of genes, RNAs, and proteins, and the `Syn3A_updated.xml` SBML file contains the metabolic reactions (reactants, stoichiometries). `The intitial_concentration.xlsx` file contains the initial count/concentrations of proteins and metabolite while `kinetic params.xlsx` contains the kinetic parameters of the GIP and metabolic reactions.


<details>
<summary><strong>Click to EXPAND: Breakdown of Input Files</strong></summary>

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

</details>

Now, we will focus on gene **JCVISYN3A_0011** to show how the different input files function in our simulation.

The genome information of JCVI-syn3A is stored as standard Genbank file on NCBI with [ACCESSION Number CP016816.2](https://www.ncbi.nlm.nih.gov/nuccore/CP016816.2/). The Genbank file holds much more information than the pure nucleotide sequence. 

Go to Genbank file under `input_data` folder and open with text editor. Search **JVCISYN3A_0011** and you wil see the following figure. **JCVISYN3A_0011** is the LocusTag of this certain gene in organism JCVI-syn3A, and *JCVISYN3A_* serving as a unique idetifier. 0011 is the locusNum that will be heavily used in the modelling to distinguish unique genes, RNAs, and proteins from each other. The start and end index of **JCVISYN3A_0011** is 15153 and 16799, respectively and the nucleotide sequence is at the end of the Genbank file. **JCVISYN3A_0011** is protein-coding gene and the protein is "Nucleoside Transporter ABC substrate-binding protein" with amino acid sequence shown. This protein is one subunit of ribonucleoside ATP-binding cassette (ABC) transporters that is assumed to import all nucleosides for Syn3A.

<p align="center">
  <img src="../figs/figs_WCM/jcvisyn3a_0011.png" width="450" alt="mpirun to launch parallel simulation"> <br>
  <b>Figure 2. Entry JCVISYN3A_0011 in Syn3A's Genbank file</b>
</p>

SBML (System Biology Markup Language) format is widely used in storing computational models of biological processes. The SBML file here contains the whole metabolic system of Syn3A, including compartments (cellular and extracellular), species, reactions, gene product (enzymatic or transporter proteins associated with reactions), and objective function (for Flux Balance Analysis). 
It's worth noting that the kinetic constants of 175 reactions are in \textit{kinetic\_params.xlsx}.

Nucleoside Transporter ABC substrate-binding protein encode by gene **JCVISYN3A_0011** functions in reaction DSGNabc, irreversible transport of deoxyguanosine in nucleotide metabolism. Figure \ref{fig:sbml} shows reaction DSGNabc in SBML file. We define \textit{reversible} as \textit{false} since DSGNabc is a irreversible transport reaction. There are three reactants and four products. The geneProductAssociation rule is *and*, meaning the four subunit proteins need to function cooperatively.

<p align="center">
  <img src="../figs/figs_WCM/sbml.png" width="450" alt="SBML entry of DSGNabc reaction"> <br>
  <b>Figure 3. Entry Reaction DSGNabc in SBML file</b>
</p>

The kinetic constants for this DSGNabc reaction in `kinetic\_params.xlsx` is shown as following.

<p align="center">
  <img src="../figs/figs_WCM/dsgnabc_kinetic.png" width="450" alt="SBML entry of DSGNabc reaction"> <br>
  <b>Figure 4. Kinetic Constants for Reaction DSGNabc</b>
</p>

## 2. Minimal Genome and Genetic Information Processes

### Genetically Minimized Bacterium, JCVI-Syn3A

Minimal cell, JCVI-syn3A is a synthetic bacterium based on mother organism Mycoplasma mycoides capri published by J. Craig Venter Institute in 2016 \cite{hutchison_design_2016}. JCVI-syn3A has a doubling time of about 2 hours and consistently forms spherical cells of approximately 200 nm in radius.

JCVI-syn3A's genome is minimal in all living organism, 543 kbp long and containing 452 genes code for proteins and 38 genes for RNAs. The ratio of unclear protein-coding genes is 90 out of 452, which is also smallest compared to other well-studied organism, including yeast and E. coli \cite{breuer_essential_2019}.

<p align="center">
  <img src="../figs/figs_WCM/syn3A_genome.png" width="450" alt="Minimized Genome"> <br>
  <b>Figure 5. Left: Protein Coding Genes of JCVI-syn3A. Only less than 90 genes are unclear.  <br>
  Right: JCVI-syn3A with genetic information processes and metabolism shown. </b>
</p>

### Genetic Information Processes

Genetic information processes connect the blueprint genes to functional proteins. In our current whole cell model of minimal cell, we mainly consider seven types of processes listed in table.

Replication copies the genetic information with the regulation of replication initiation. Transcription copies sequential information from DNA to RNA. There are three types of RNA, including mRNA, rRNA and tRNA that are tightly connected by translation. Translation takes place on ribosomes, where mRNA is read and an amino acid chain is generated according to the sequence of mRNA. rRNA, together with other ribosome proteins make up ribosomes in ribosomal biogenesis \cite{earnest_ribosome_2016}. tRNA is the carriers and identifier of amino acid to the ribosome. 

<p align="center">
  <img src="../figs/figs_WCM/gip_table.png" width="450" alt="GIP reactions"> <br>
  <b>Figure 6. Left: Genetic Information Flow. <br>
  Right: Functions of GIP reactions. </b> 
</p>


The translation and degradation reactions in current whole cell model are more precise compared to the simplest case in Tutorial 2 as for each reaction, mRNA first needs to bind with a complex machinery, degradosome or ribosome. As the busiest species in genetic information processes, mRNA can also be degraded by binding with degradosome. This competition of mRNA to bind with ribosome or degradosome is a important pathway to regulate genetic information processes that will shown in the following result.


Replication initiation is modeled by the binding of multi-domain protein DnaA (encoded by JCVISYN3A\_0001) and replisome with certain region called oriC on the chromosome \cite{thornburg_kinetic_2019}. 

<p align="center">
  <img src="../figs/figs_WCM/oric_dnaA.png" width="450" alt="DnaA with Ori"> <br>
  <b>Figure 7. Left: Left: oriC region. 9 nucleotide signature binding with DnaA domain IV shown in yellow and red, 3 nuclotide AT rich region binding with DnaA domain III shown in grey. <br>
  Right: PDB structure of DnaA domain IV and domain III binding with chromosome a): \cite{fujikawa_structural_2003} b): \cite{duderstadt_dna_2011} </b> 
</p>

First stage of initiation is the binding of DnaA's Domain IV with nine-nucleotide signatures on double-strand DNA(three of them, shown in red and yellow in Figure \ref{oric_dnaA}). This binding will open a pocket for DnaA's Domain III to bind with AT rich region on single-strand DNA following the nine-nucleotide signatures to build a filament of DnaA on single-strand DNA. The final step is the binding of replisome with the chromosome after the filament grows above 15 DnaA with 30 DnaA maximum due to the length of AT rich region is 90. In current model, the assembly of replisome is not explicitly included.

<p align="center">
  <img src="../figs/figs_WCM/initiation.png" width="450" alt="DnaA with Ori"> <br>
  <b>Figure 8. Three stage DNA replication initiation </b> 
</p>

Replication occurs gene by gene when replisomes move along the circular chromosome from oriC to ter in two directions. Since the locations of genes are fixed, the replication happens in order. So our replication model is ordered and bidirectional.

<p align="center">
  <img src="../figs/figs_WCM/rep.png" width="450" alt="DnaA with Ori"> <br>
  <b>Figure 9. Ordered and Bidirectional Replication </b> 
</p>

The equivalent reaction to duplicate certain gene is in Table \ref{tab:CMErxns} The rate to replicate one gene in Table \ref{tab:CMErates} is determined by a polymerization rate form proposed in \cite{hofmeyr_generic_2013}.  

**Table 1: Reactions and Rates for Replication, Transcription, Translation and Degradation**

| **Processes**  | **Reactions**                                                                 | **Kinetic Constant**                      |
|----------------|-------------------------------------------------------------------------------|-------------------------------------------|
| Replication     | $G_{\text{locusNum}} \rightarrow 2G_{\text{elongation}}$                      | $k_{\text{replication}}^{\text{locusNum}}$ |
| Transcription   | $RNAP + G_{\text{locusNum}} \rightarrow RNAP:G_{\text{locusNum}}$            | $k_{\text{trsc}}^{\text{binding}}$         |
|                | $RNAP:G_{\text{locusNum}} \rightarrow R_{\text{locusNum}} + RNAP + G_{\text{locusNum}}$ | $s_{\text{trsc}}^{\text{locusNum}} k_{\text{trsc}}^{\text{elongation}}$ |
| Translation     | $Ribosome + R_{\text{locusNum}} \rightarrow Ribosome:R_{\text{locusNum}}$     | $k_{\text{trans}}^{\text{binding}}$        |
|                | $Ribosome:R_{\text{locusNum}} \rightarrow P_{\text{locusNum}} + Ribosome + R_{\text{locusNum}}$ | $k_{\text{trans}}^{\text{elongation}}$ |
| Degradation     | $Degradosome + R_{\text{locusNum}} \rightarrow Degradosome:R_{\text{locusNum}}$ | $k_{\text{degra}}^{\text{binding}}$       |
|                | $Degradosome:R_{\text{locusNum}} \rightarrow NMPs + Degradosome$              | $k_{\text{degra}}^{\text{depoly}}$        |

Transcription, translation and degradation are all depicted as a two-step binding and (de)polymerizatoin reactions where the polymerization in transcription and translation shares the same rate form as that in replication. The rate for degradation from a mRNA to its monomers is calculate by divide the monomer depletion rate over the length of mRNA.

**Table 2: Rate Form for Replication, Transcription, Translation and Degradation**

| **Processes**   | **(De)polymerization Constants**                  | **Formulation**                                                                                                                                  |
|-----------------|---------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| Replication     | $k_{\text{elongation}}^{\text{replication}}$      | $\dfrac{k_{\text{replication}}^{\text{cat}}}{\dfrac{K_{D1} K_{D2}}{[dNTP_1][dNTP_2]} + \sum_i \dfrac{K_{Di}}{[dNTP_i]} + L_{\text{DNA}} - 1}$     |
| Transcription   | $k_{\text{elongation}}^{\text{trsc}}$             | $\dfrac{k_{\text{transcription}}^{\text{cat}}}{\dfrac{K_{D1} K_{D2}}{[NTP_1][NTP_2]} + \sum_i \dfrac{K_{Di}}{[NTP_i]} + L_{\text{RNA}} - 1}$     |
| Translation     | $k_{\text{elongation}}^{\text{trans}}$            | $\dfrac{k_{\text{translation}}^{\text{cat}}}{\dfrac{K_{D1} K_{D2}}{[\text{tRNA:aa}_1][\text{tRNA:aa}_2]} + \sum_i \dfrac{K_{Di}}{[\text{tRNA:aa}_i]} + L_{\text{protein}} - 1}$ |
| Degradation     | $k_{\text{depoly}}^{\text{degra}}$                | $\dfrac{k_{\text{cat}}^{\text{degra}}}{L_{\text{mRNA}}}$                                                                                        |

## 3. Essential Metabolism and Rate Law

### Essential Metabolism

As mentioned, replication, transcription, and translation all require monomers (deoxyribonucleotide (dNTPs), ribonucleotide (NTPs), and amino acids (AAs)) for the polymerizations reactions. 

All the supply of monomers comes from metabolism of JCVI-syn3A. The (d)NTPs are generated in nucleotide metabolism that transport extracellular nucleoside and convert nucleoside to (d)NTPs used in the DNA and RNA synthesis. The other building blocks, such as phosphoribosyl pyrophosphate (prpp), phosphoenolpyruvate (pep), and 1,3-Bisphosphoglyceric acid (1,3dpg) along these pathways are the products in glycolysis in central metabolism \ref{NTP_supply}.

<p align="center">
  <img src="../figs/figs_WCM/NTP_supply.png" width="450" alt="GTP synthesis"> <br>
  <b>Figure 10. Supply of GTP and dGTP for DNA and RNA synthesis in nucleotide metabolism. The metabolites in the bracket come from central metabolism </b> 
</p>

The whole metabolism network can be separated to five connected sub-networks, including central, nucleotide, lipid, cofactor and amino acid. Glucose 6-phosphate (g6p) the immediate downstream of glucose in glycolysis connects central and lipid metabolism, while Adenosine triphosphate (ATP) that generated in glycolysis provides energy for reactions in all five sub-networks. 

<p align="center">
  <img src="../figs/figs_WCM/metabolism.png" width="450" alt="GTP synthesis"> <br>
  <b>Figure 11. Whole Metabolism of JCVI-syn3A with five subsystems: central, nucleotide, lipid, cofactor and amino acid. </b> 
</p>

In our current whole cell modeling, 197 metabolites are included, 148 cytoplasmic and 49 extracellular. 175 reactions connect the metabolites. Shown in the network of central metabolism, the orange nodes are the metabolite IDs. We use suffix '\_c' to denote cytoplasmic and '\_e' for extracellular. The blue arrows are reactions connecting different metabolites, with names in dark blue and related proteins' locusNums under the names. 

ATP as the major energetic molecules energize cellular processes, without which the cell will die. In JCVI-syn3A, glycolysis is the only way to generate ATP. Phosphoenolpyruvate (pep) is an intermediate in both the upstream and downstream in glycolysis.

<p align="center">
  <img src="../figs/figs_WCM/central_pep.png" width="450" alt="GTP synthesis"> <br>
  <b>Figure 12. Central Metabolism in JCVI-syn3A. Phosphoenolpyruvate (pep) is both the upstream and downstream of glycolysis. 1,3-Bisphosphoglyceric acid (1,3dpg), Phosphoribosyl pyrophosphate (prpp) and pep are in nucleotide metabolism. </b> 
</p>


### Convenience Rate Law

Most reactions in the metabolism requires certain proteins produced in genetic information processes as enzymes or transporters. Convenience rate law and random binding model are used to depict the kinetics of these reactions \cite{liebermeister_bringing_2006}.

The random binding model assumes the substrates bind to the enzyme/transporter in arbitrary order and are converted into the products, which then dissociate from the enzyme in arbitrary order. Convenience rate law further assumes the conversion step is the rate-determining step and the binding reactions are in quasi-equilibrium. For reaction

```math
\ce{A + X <=>[k_f][k_r] B + Y}
```

, the rate law is

<p align="center">
  <img src="../figs/figs_WCM/convenience.png" width="450" alt="GTP synthesis"> <br>
  <b>Figure 13. Mechanism of reaction A+X to B+Y in random binding model. [E] is the concentration of enzyme, [A] concentration of molecule A </b> 
</p>

The metabolic network is also stored in standard Systems Biology Markup Language (SBML) file. In the given SBML file of JCVI-syn3A, there are two compartments, cellular and extracellular with metabolites and reactions. However, SBML file itself does not store kinetic parameters and we choose to store them in a Excel file with multiple sheets for the simulation.

We can look at one reaction entry in the SBML file for more insights. Go to the SBML file and search DGSNabc. DGSNabc, the irreversible transport of deoxyguanosine (one of RNS ABC import system) in nuclotide metabolism for example, has three reactants and four products with four proteins involved. The geneProductAssociation is AND, meaning four proteins are all needed to perform this reactions. Biologically, four proteins are four sub-units of nucleoside ABC transporter.


## 4. Hybrid CME-ODE Algorithm

### Step-wise Communication Between CME and ODE

The discreteness and stochasticity of chemical kinetics play a role when the number of reactants is significantly low. This makes it necessary to use stochastic chemical master equation (CME) for genetic information processes (GIP) where the copy numbers of species are low, shown as follows to sample the variation. In contrast, ordinary differential equations (ODE) is sufficient to depict the kinetics of large numbers of small metabolites in metabolism.

$\frac{dP(\mathbf{x},t)}{dt}=\sum_{r}^{R} [-a_r({{\mathbf{x}}}) P({{\mathbf{x}}},t) + a_r({{\mathbf{x}}}_\nu-\mathbf{S_r}) P({{\mathbf{x}}}-\mathbf{S_r},t)]$

To simulate the **co-evolution** of GIP and metabolism, the communication needs to be performed to describe the interactions between these two subsystems. We then proposed to use **[Hybrid CMEODE simulation](https://ietresearch.onlinelibrary.wiley.com/doi/10.1049/iet-syb.2017.0070)**.  We first discretize the entire simulation length into piecewise communication time steps (hook intervals, $t_H$). During each communication time step, 

**(a)** A CME simulation of length $t_H$ is performed to describe the kinetics in GIP. 

**(b)** The communication from CME to ODE by passing the protein counts, the consumption of monomers (dNTP, NTP, and amino acid charged tRNAs), and the recycling of nucleoside monophosphates (NMP) in the degradation of mRNA and recycling of amino acids in the degradation of membrane protein. 

**(c)** Then a $t_H$ length ODE simulation is performed with the updated concentrations of proteins and metabolites. 

**(d)** The impacts of metabolism on GIP are two-fold: the abundance of metabolites explicitly in GIP and the concentrations of monomers that affect the rates in the polymerization of gene, RNA, and protein.

<p align="center">
  <img src="../figs/figs_WCM/communication.png" alt="GIP_Cplx" width="400">
</p>

### Program Flowchart

The CME simulation is executed using Lattice Microbes (LM) with direct Gillespie algorithm. We employ the `hookSimulation` function to interrupt the CME timeline and enable communication with the ODE solver. For the ODE simulation, we use the **[odecell](https://github.com/Luthey-Schulten-Lab/odecell)** software developed by the Luthey-Schulten Lab, which maps metabolic reactions to ordinary differential equations and specifies the corresponding kinetic parameters. The resulting ODE system is solved using the *lsoda* algorithm from the SciPy library.

One extra thing to notice is that the CME rates are also updated per second after the metabolism simulation in ODE to account for the possible shortage of monomers that will decrease the rate in genetic information process.

<p align="center">
  <img src="../figs/figs_WCM/Flowchart_CMEODE_WCM_Poster.png" alt="GIP_Cplx" width="400">
</p>

## 5. Analysis

### Run ~~pkl.ipynb~~ then plotting.ipynb on Jupyter Notebook Webpage
+ **First**: Navigate to `./LM/CME/WholeCellModel/analysis` in Jupyter Notebook Webpage.
  
+ **Second**: Open *pkl.ipynb* and **Run All** to serialize the CSV files into one single python pickle file.  
  10 CSV files are generated beforehand for the analysis.
  
+ **Thrid**: Open *plotting.ipynb* and **Run All** to generate the plots.  
  Compare the plots with figures in *Summer_School_CME_2024.pdf*.

### Cell growth

### Doubling of Proteome

### Traces of important metabolite

## 6. Discussion

### DNA initiation

### Stochastic translation of protein

### Translation per mRNA

### Slowdown of GIP