#!/usr/bin/env python
# coding: utf-8

# # Tutorial 1.2 - Stochastic Genetic Information Process
# 
# Here we examine a CME model of stochastic Genetic Information Process.
# 
# In this model, we include the transcription and translation of gene and mRNA together with degradation of both mRNA and protein.
# 
# The model presented here can be found in the article: [Analytical distributions for stochastic gene expression](https://www.pnas.org/doi/full/10.1073/pnas.0803850105).
# 

# In[ ]:


# Import Standard Python Libraries
import os
import numpy as np
import matplotlib.pyplot as plt

# Import pyLM Libraries
from pyLM import *
from pyLM.units import *
from pySTDLM import *
from pySTDLM.PostProcessing import *

# Enable plotting inline in the Jupyter notebook
get_ipython().run_line_magic('matplotlib', 'inline')


# ## Constants
# 
# Rates of reactions come from the [cell](https://www.cell.com/cell/fulltext/S0092-8674(21)01488-4?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421014884%3Fshowall%3Dtrue) paper's Whole Cell Model of DnaA coding gene (G_0001) at initial conditions.
# 
# Degradation rate of protein is calculated based on 25 hours' half life in [Maier et al, 2011](https://www.embopress.org/doi/full/10.1038/msb.2011.38).

# In[ ]:


# Constants
k_transcription  = 6.41e-4       # Transcription, s^-1
k_degra_mRNA = 2.59e-3     # degradation of mRNA, s^-1
k_translation = 7.2e-2        # translation, s^-1
k_degra_ptn = 7.70e-6      # degradation of protein, s^-1


# ## Define CME simulation

# We begin by creating a [CMESimulation](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) "object" that we call ```sim```. This object will include the definition of the whole stochastic simulation.

# In[ ]:


# Create our CME simulation object
sim = CME.CMESimulation(name='Gene Expression')


# Next we define the chemical species with simulation. First. we specify the names of the chemical species.  Then we register these species with the simulation.  The [```defineSpecies()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) function can be called multiple times and will add any new names to the list of species.

# In[ ]:


# Define our chemical species
species = ['gene', 'mRNA', 'ptn']
sim.defineSpecies(species)


# Here we add reactions to the simulation. We use the [```addReaction()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) function that is a member of the ```CMESimulation``` object. We add a bimolecular association reaction and a unimolecular dissociation reaction. When more than one reactant is involved, the list of reactant names should be passed as a tuple as can be seen in the reactant of the association reaction. The rates in this command must be in units of molecules and seconds, for instance units of ```/molecule/sec``` for the association reaction.

# In[ ]:


# Add reactions to the simulation

sim.addReaction(reactant='gene', product=('gene','mRNA'), rate=k_transcription)
sim.addReaction(reactant='mRNA', product='', rate=k_degra_mRNA)
sim.addReaction(reactant='mRNA', product=('mRNA','ptn'), rate=k_translation)
sim.addReaction(reactant='ptn', product='', rate=k_degra_ptn)


# Next, we add the initial particle counts to the simulation using the [```addParticles()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) function.

# In[ ]:


# Set our initial species counts

sim.addParticles(species='gene', count=1)
sim.addParticles(species='mRNA', count=1)
sim.addParticles(species='ptn', count=0)


# Finally, we define the simulation execution parameters. We have the simulation run for 6300 seconds of real time to cover the entire cell cyle.
# 
# The traces are recorded per 1 second.
# 
# Then we name the simulation output file and save the simulation definition to it.

# In[ ]:


# Simulation time is 6300 (entire cell life cycle of Minimal Cell).
writeInterval = 1
simtime = 6300

sim.setWriteInterval(writeInterval)
sim.setSimulationTime(simtime)

filename = "./T2.1-GeneticInformationProcess.lm"

os.system("rm -rf %s"%(filename)) # Remove previous LM file 

sim.save(filename)


# In[ ]:


# Print out the information of the system
sim


# ## Run Simulation

# In[ ]:


# Run multiple replicates using the Gillespie solver
reps = 100

sim.run(filename=filename, method="lm::cme::GillespieDSolver", replicates=reps)


# ## Post-Processing

# Create Picture Folder

# In[ ]:


plotfolder = './plots_GeneticInformationProcess/'

if not os.path.exists(plotfolder):
    os.mkdir(plotfolder)


# Using function [```plotTraceFromFile```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pySTDLM.PostProcessing.html#pySTDLM.PostProcessing.plotTraceFromFile) to plot the trace of mRNA or protein in the specified simulation replicate.

# In[ ]:


rep = 3

picturepath = plotfolder + 'mRNA_replicate{0}.png'.format(rep)

PostProcessing.plotTraceFromFile(filename, species=['mRNA'], replicate=rep, outfile=picturepath )

picturepath = plotfolder + 'ptn_replicate{0}.png'.format(rep)

PostProcessing.plotTraceFromFile(filename, species=['ptn'], replicate=rep, outfile=picturepath )


# Using [```plotAvgVarFromFile()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pySTDLM.PostProcessing.html#module-pySTDLM.PostProcessing) built-in function in PostProcessing to plot the average and variance of all three species. 

# In[ ]:


species = ['gene', 'mRNA', 'ptn']

for i_specie, specie in enumerate(species):
    
    if specie != 'gene':
        
        picturepath = plotfolder + 'AvgVar_{0}_{1}replicates.png'.format(specie, reps)
        
        PostProcessing.plotAvgVarFromFile(filename=filename, species=[specie], outfile=picturepath)


# #### Using h5py module to serialize traces in LM file to a 3D Numpy Array with dimesions *(reps, species, time)*.

# In[ ]:


import h5py

fileHandle = PostProcessing.openLMFile(filename) # Create h5py file handle

timestep = PostProcessing.getTimesteps(fileHandle) # use PostProcessing to get the timesteps of the simulation

traces = np.zeros((reps,len(sim.particleMap),len(timestep)))

def get_sim_data(filename):

    f=h5py.File(filename, "r")

    for r in range(reps):
        
        traces[r]=f['Simulations'][str(r+1).zfill(7)]['SpeciesCounts'][()].transpose()
        
    f.close()

    return traces

traces = get_sim_data(filename)

print('The size of the 3D trajectories is {0} with dimensions replicates, species, and time.'.format(np.shape(traces)))


# #### Plot the distribution of Protein at the end of the simulation

# In[ ]:


import plot_hist# User defined function to draw histogram

# Protein Distribution

ptns_end = traces[:,2,-1] # Slice Numpy array to get the counts of ptn at the end of the whole cell cycle

fig_path = plotfolder + 'Distribution_Ptns_at_{0}_seconds_{1}replicates.png'.format(simtime, reps)

xlabel = 'Counts of Protein [#]'

title = 'Distribution of Ptn Counts at {0} seconds {1} replicates'.format(simtime, reps)
    
plot_hist.plot_histogram(data=ptns_end, figure_path=fig_path, bins=20, xlabel=xlabel, title=title)

# mRNA Distribution
###################################################################################
mRNAs_end = traces[:,1,-1] # Slice Numpy array to get the counts of ptn at the end of the whole cell cycle

fig_path = plotfolder + 'Distribution_mRNA_at_{0}_seconds_{1}replicates.png'.format(simtime, reps)

xlabel = 'Counts of mRNA [#]'

title = 'Distribution of mRNA Counts at {0} seconds {1} replicates'.format(simtime, reps)

plot_hist.plot_histogram(data=mRNAs_end, figure_path=fig_path, bins=5, xlabel=xlabel, title=title)


# ## Questions 2.1
# 
# 1. Do mRNA and protein reach steady-state during the 6300 seconds' simulation? How can you tell this from the plots?
# 
# 2. The initial count of protein P\_0001/DnaA from experimental proteomics data is 148. Compare the mean count of protein at the end of the cell cycle to this experimental count. Does the simulation roughly generate 148 proteins during the entire cell cycle? And why this is important?
# 
# 3. \(Challenge:\) What should be the distribution of protein counts in different replicates? You can find an analytical solution in [Swain's paper](https://www.pnas.org/doi/full/10.1073/pnas.0803850105) in 2008. The distribution of steady state protein distribution is ploted in the following block. Why the two distributions are so different? Does protein count reach steady state during the simulation?

# ### Negative binomial distribution of protein counts
# 
# According to Swain's paper, the steady state distribution of protein count should be negative binomial distribution. The following code compares the simulation result with analytical one.

# In[ ]:


a = k_transcription/k_degra_ptn # Number of mRNA in a protein's lifetime
b = k_translation/k_degra_mRNA # b/(1+b) is the probability that one mRNA is translated rather than degradated

print("Number of mRNA in a protein's lifetime is {0}".format(a))
print('THe probability that one mRNA is translated rather than degradated is {0}'.format(b/(1+b)))


# In[ ]:


import numpy as np
from scipy.stats import nbinom
import matplotlib.pyplot as plt

# Enable plotting inline in the Jupyter notebook
get_ipython().run_line_magic('matplotlib', 'inline')

fig_size = [87,87/1.618]

fig, ax = plt.subplots(1, 1)

n, p = int(np.floor(a)), np.round(1/(1+b),2) 

x = np.arange(nbinom.ppf(0.01, n, p),
              nbinom.ppf(0.99, n, p))

plt.title(label = 'Analytical Distribution of Steay State Protein Count')

ax.plot(x, nbinom.pmf(x, n, p),color='c', label='Negative Binomial Distribution')

x_mean = np.sum(x*nbinom.pmf(x, n, p)) # Compute the mean by summing up the product of probability and x

ax.vlines(x_mean, 0, np.max(nbinom.pmf(x,n,p)), linestyles='dashed',colors='b', lw=2, alpha=1, label='Mean: {0:.0f}'.format(x_mean))

ax.legend(loc='upper left', fontsize='large',frameon=False)

plt.savefig(plotfolder + 'Distribution_Steady_State_ptn.png')

plt.show()

