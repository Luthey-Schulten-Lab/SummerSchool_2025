# Stochastic Genetic Information Processs in CME

Open [ChatGPT](https://chatgpt.com/) and ask this question:

**Give an example on large differences in rate constants that lead to significant fluctuations**

You can ask ChatGPT to give you more examples. Please also check this answer: [stochastic genetic information processes](https://chatgpt.com/c/91d87e61-bae5-4b89-b078-d1d4cfa44274), which is this tutorial.

## Classic Genetic Information Processes Reactions

With 3 species and 4 reactions, this classic and simplest genetic information process (GIP) starts from the transcription of gene to mRNA. mRNA can be translated to protein or degraded to its monomers. Protein can also be degraded. The reactions and rate constant is shown as the followings.

<p align="center">
  <img src="../figs/figs_GIP/GIP_withCMEs.png" width="600" alt="Simple GIP model">  <br>
  <b>Figure 1. Genetic information processing model and its chemical master equation, where <i>m, n</i> are the number of mRNAs and proteins </b>
</p>

The rate constants are for DnaA Coding Gene (G\_0001) of minimal cell. The first three rate constants are calculate based on the initial concentrations of nucleotides and amino acids charged tRNA in the Cell paper[^thornburg_cell]. We also fix the gene copy number to 1 and assume initial count of mRNA to be 1. The initial count of protein is 0. The degradation rate of protein is estimated based half-life 25 hours[^thornburg_kinetic].

**Table 1. Four reactions with their rate constants**

| **Names**              | **Reaction**                          | **Rate Constant (s<sup>-1</sup>)**                              | **Propensity (s<sup>-1</sup>)**                              |
|------------------------|----------------------------------------|------------------------------------------------------|---------------------------------------------------|
| Transcription          | Gene → mRNA                            | *k*<sub>transcription</sub> = 6.41×10<sup>-4</sup>   | *k*<sub>transcription</sub>                       |
| Degradation of mRNA    | mRNA → ∅                               | *k*<sub>deg,m</sub> = 2.59×10<sup>-3</sup>           | *k*<sub>deg,m</sub> · *N*<sub>mRNA</sub>          |
| Translation            | mRNA → mRNA + Protein                  | *k*<sub>translation</sub> = 7.20×10<sup>-2</sup>     | *k*<sub>translation</sub> · *N*<sub>mRNA</sub>    |
| Degradation of Protein | Protein → ∅                            | *k*<sub>deg,p</sub> = 7.70×10<sup>-6</sup>           | *k*<sub>deg,p</sub> · *N*<sub>ptn</sub>           |

## Run the Jupyter Notebook

Now go to file on your Jupyter Notebook webpage `Tut2.1-GeneticInformationProcess.ipynb` to simulate this toy GIP model. The default simulation length `simtime` is 6300 seconds, the entire cell cycle of minimal cell.We simulate 10 independent cell `reps` and write out the trajectories at `writeInterval` of 1 second.

## Stochastic Protein Synthesis

We first look at the average and span of mRNA and protein abundances among the cell population of 10 replicates. The population averaged mRNA abundace fluctuates below 1 along the cell cycle, and the protein counts increase accumulatively from the translation processes (protein degradation was minor compared to translation under this set of kinetic parameters).

Protein synthesis takes place in each single cells. We can plot the traces of mRNA and protein in different single replicates to see the stair-stepping trace of mRNA and the burst of protein. You will an increase/burst in protein when there are mRNAs, and the halting or even decrease of protein with no mRNA. You are encouraged to compare the pattern shown here to your own plots.

<p align="center">
  <img src="../figs/plots_GIP/GIP_mRNA_Protein_10Replicates.png" width="300" alt="mRNA Protein 10 reps"> <img src="../figs/plots_GIP/GIP_mRNA_Protein_Cell1.png" width="300" alt="CME replicate 1"> <br>
  <b>Figure 2. Left: Population average (solid line) and full span (shaded area) of mRNA (red) and protein(blue) abundances in 10 cell replicates <br> 
  Right: star-stepping trace of mRNA and the burst-like protein synthesis in one single cell replicate</b>
</p>

## Discussion

### 1. Stedy-state

Do mRNA and protein reach steady-state during the 6300 seconds' simulation? How can you tell this from the plots? If the fluctuation is large, try to increase the replicates numbers `reps` from 10 to 100.

### 2. Doubling the initial abundace of protein for cell division
The initial count of protein P\_0001/DnaA from experimental proteomics data is 148. In the following histogram, the population average of DnaA at the end of the cell cycle is 270. Compare the mean count of protein at the end of the cell cycle to this experimental count. Does the simulation roughly generate 148 proteins during the entire cell cycle? And why this is important? Please consider cell division.

<p align="center">
  <img src="../figs/plots_GIP/GIP_Proteins_CycleEnd_100replicates.png" width="450" alt="ODE result"> <br>
  <b>Figure 3. Distribution of protein abundances among 100 cell replicates at the of the cell cycle</b>
</p>

## References:
[^thornburg_cell]: Thornburg, Z. R., Bianchi, D. M., Brier, T. A., Gilbert, B. R., Earnest, T. M., Melo, M. C., Safronova, N., Sáenz, J. P., Cook, A. T., Wise, K. S., Hutchison, C. A., Smith, H. O., Glass, J. I., & Luthey-Schulten, Z. (2022). Fundamental behaviors emerge from simulations of a living minimal cell. Cell, 185(2), 345-360.e28. https://doi.org/10.1016/j.cell.2021.12.025

[^thornburg_kinetic]: Thornburg, Z. R., Melo, M. C. R., Bianchi, D., Brier, T. A., Crotty, C., Breuer, M., Smith, H. O., Hutchison, C. A., Glass, J. I., & Luthey-Schulten, Z. (2019). Kinetic modeling of the genetic information processes in a minimal cell. Frontiers in Molecular Biosciences, 6. https://doi.org/10.3389/fmolb.2019.00130
