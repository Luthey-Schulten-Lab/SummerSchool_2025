#!/usr/bin/env python
# coding: utf-8

# # Tutorial 1.2 - CME Solution of a Bimolecular Reaction
# 
# ### Here we examine a stochastic version of Tutorial 1.1 using Chemical Master Equation.
# 
# We use pyLM to construct the CME system.
# 
# In order to use pyLM, we need to import several libraries.  The first is [```pyLM proper (pyLM)```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/API.html).  The second [```pyLM.units```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.units.html#module-pyLM.units) is a library to convert number with given units, such as ```nm()```, ```micron()```, ```ms()```, ```microsecond()```, etc.  Finally, we import the [```pySTDLM```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pySTDLM.html#module-pySTDLM) standard library, which contains functionality used in biological simulations and post-processing commands including plotting. [```pySTDLM.PostProcessing```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pySTDLM.PostProcessing.html#module-pySTDLM.PostProcessing) module contains functions to get traces and plotting.

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


# In[ ]:


# Constants
V  = 1.0e-15       # L
NA = 6.022e23      # molecules/mole

kf_ODE = 1.07e6; # /M/s
kr_ODE = 0.351 # /s

kf = 1.07e6/(NA*V) # convert from 1.07e5 /M/s to /counts/s
kr = 0.351         # /s


# ## Define a CME Simulation
# 
# See the [pyLM Documentation](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/API.html)
# 
# We begin by creating a [CMESimulation](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) "object" that we call ```sim```. This object will include the definition of the whole stochastic simulation.

# In[ ]:


# Create our CME simulation object
sim = CME.CMESimulation(name='Bimolecule Reaction')


# Next we define the chemical species with simulation. First. we specify the names of the chemical species.  Then we register these species with the simulation.  The [```defineSpecies()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) function can be called multiple times and will add any new names to the list of species.
# 

# In[ ]:


# Define our chemical species
species = ['A', 'B', 'C']
sim.defineSpecies(species)


# Here we add reactions to the simulation. We use the [```addReaction()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) function that is a member of the ```CMESimulation``` object. We add a bimolecular association reaction and a unimolecular dissociation reaction. When more than one reactant is involved, the list of reactant names should be passed as a tuple as can be seen in the reactant of the association reaction. The rates in this command must be in units of molecules and seconds, for instance units of ```/molecule/sec``` for the association reaction.
# 

# In[ ]:


# Add reactions to the simulation
sim.addReaction(reactant=('A','B'), product='C', rate=kf)
sim.addReaction(reactant='C', product=('A','B'), rate=kr)


# Next, we add the initial particle counts to the simulation using the [```addParticles()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pyLM.CME.html#module-pyLM.CME) function.

# In[ ]:


# Set our initial species counts
sim.addParticles(species='A', count=50)
sim.addParticles(species='B', count=50)
sim.addParticles(species='C', count=0)


# Finally, we define the simulation execution parameters. We have the simulation run for 30 seconds of real time saving results every 30 microseconds for a total of 1 million times. <br/>
# 
# Then we name the simulation output file and save the simulation definition to it.

# In[ ]:


# Define simulation parameters: run for 30 seconds, saving data every 30 ms
sim.setWriteInterval(microsecond(30))
sim.setSimulationTime(30)
filename = "./T1.2-bimol.lm"
os.system("rm -rf %s"%(filename))
sim.save(filename)


# In[ ]:


# Print simulation parameters to the notebook
sim


# ## Run Simulation with Different Replicates Number
# 
# Next we run the simulation. To do this, we specify which file has the problem specification (as saved two cells up). Lattice Microbes supports several different solvers for CME simulations; here we use the most common algorithm called the Gillespie's algorithm (aka the Stochastic Simulation Algorithm). You can learn more from Gillespie's [review paper](https://www.annualreviews.org/content/journals/10.1146/annurev.physchem.58.032806.104637) and [wikipedia page](https://en.wikipedia.org/wiki/Gillespie_algorithm).
# 
# Because the CME represents a stochastic process, each instance of a simulation will have a slightly different trajectory.  Generally, we run many "replicate" simulations with the same or nearly the same starting conditions and then compute aggregate statistics over all of them.

# In[ ]:


# Run reps replicates using the Gillespie solver
reps = 10

sim.run(filename=filename, method="lm::cme::GillespieDSolver", replicates=reps)


# ## Post-Process Simulation
# 
# See the [pySTDLM Documentation](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pySTDLM.PostProcessing.html#module-pySTDLM.PostProcessing)
# 
# Post-processing generally begins by getting a handle to the file. This is accomplished by passing the filename to the function [```openLMFile()```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pySTDLM.PostProcessing.html#pySTDLM.PostProcessing.openLMFile). This function does some error checking to make sure the file is generated by LM. <br/>
# 
# The function [```plotTrace```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pySTDLM.PostProcessing.html#pySTDLM.PostProcessing.plotTrace) will plot a list of species from the specified simulation replicate.
# 
# Finally, we close the LM File. It is very important that if you open an LM file with the function ```openLMFile``` that it be closed at the end of your post-processing with [```closeLMFile```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pySTDLM.PostProcessing.html#pySTDLM.PostProcessing.closeLMFile).   The function takes the filehandle that is returned by ```openLMFile``` as an argument. Also, not that once this function is called, any function that takes the filehandle as an argument will fail to work as the handle is now stale.  This is a common mistake and if you get crashing, check that you haven't prematurely closed the file.  This function is usually called last in a post-processing script.
# 
# 

# In[ ]:


plotfolder = './plots_bimolecule/'

if not os.path.exists(plotfolder):
    os.mkdir(plotfolder)


# Plot the traces of one replicate to see the fluctuation

# In[ ]:


fieHandle = PostProcessing.openLMFile(filename)

rep = 1

picturepath = plotfolder + 'bimolecule_CME_trace_replicate{0}.png'.format(rep)

PostProcessing.plotTrace(fieHandle, species=['A','C'], replicate=1, filename=picturepath)

closeLMFile(fieHandle)


# The function [```ploAvgVarFromFile```](https://luthey-schulten.chemistry.illinois.edu/software/LM2.4/_autosummary/pySTDLM.PostProcessing.html#pySTDLM.PostProcessing.plotAvgVarFromFile) will plot averages and variances of a list of species in one plot for all replicates. The plot will be saved in path ```outfile```. 

# In[ ]:


# Plot average and variance among all replicates
picturepath = plotfolder + 'bimolecule_CME_Avgs_Vars_{0}replicates.png'.format(reps)

PostProcessing.plotAvgVarFromFile(filename = filename, species = ['A','C'], outfile = picturepath)


# ## Questions 1.2
# 
# 1. Try plotting the average and variance of different replicates by changing ***reps*** from 10 to 100 even more. You need to restart the kernal to start a new CME simulation.
# 
# 2. How many replicates are required to get a smooth average? How many for a smooth variance?
# 
# 4. (Challenge) Can you derive an analytical solution for the system of equations?  Try fitting the rate constants using the results of the stochastic simulations. You may find [```scipy.optimize.curve_fit()```](http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html) useful for this.
# 
# 5. (Challenge) Does the theoretical average number from CME always be consistent/same with the count from ODE? You may find Page 420 in McQuarrie's classic paper [STOCHASTIC APPROACH TO CHEMICAL KINETICS](https://edisciplinas.usp.br/pluginfile.php/5860923/mod_resource/content/2/McQuarrie_Gillespie.pdf) helpful. 
# 
# 

# 
