"""
User defined Python Script to use matplotlib to plot the distribution of protein and mRNA count
"""
import os
import numpy as np
import h5py 
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps



mm = 1/25.4
l_w_ratio = 1/1.618


def get_sim_data(traces, reps, filename):
    """
    Description: serialize the trajectories from the output LM file into a 3D numpy array
    """


    f=h5py.File(filename, "r")

    for r in range(reps):

        traces[:,:,r]=f['Simulations'][str(r+1).zfill(7)]['SpeciesCounts'][()]

    f.close()


    return traces

def choose_colormap(N_items):
    """
    Description: Set the color map
    """

    if N_items <= 10:
        cmap = colormaps.get_cmap('tab10')
        c_space = np.arange(0,N_items,dtype=np.int32)
    else:
        cmap = colormaps.get_cmap('plasma')
        c_space_lower_lim = 0.0
        c_space_upper_lim = 0.9
        c_space = np.linspace(c_space_lower_lim,
                                c_space_upper_lim,
                                N_items,
                                dtype= np.float32)
    return cmap, c_space


def plot_traces(timestep, traces, legends, fig_path, ylabel, title):
    """
    Description: plot the time traces of different species in single replicate
    
    """
    fig_size = [87,87*l_w_ratio]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'time [s]',
                    fontsize=7,
                    labelpad=1.5)
    ax.set_ylabel(r''+ylabel.replace('_','_'),
                fontsize=7,
                labelpad=1.5)
    ax.set_title(r''+title.replace('_','_'), 
                fontsize=8,
                pad=4)
    ax.set_xlim(timestep[0] - 0.05*timestep[-1],
                timestep[-1])

    ax.set_ylim(np.amin(traces),
                1.2*np.amax(traces))
    
    tick_length = 2.0
    tick_width = 1.0
    ax.tick_params(labelsize=5,
                    length=tick_length,
                    width=tick_width,
                    direction='in',
                    left=True,
                    right=True,
                    bottom=True,
                    top=True,
                    which='major')


    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)

    n_species = np.shape(traces)[1]
    cmap, c_space = choose_colormap(n_species)

    for i_specie in range(n_species):
        trace_specie = traces[:,i_specie]
        ax.plot(timestep,trace_specie,
                                alpha=0.5,
                                linewidth=0.5,
                                zorder=-3,
                                color=cmap(c_space[i_specie]),
                                label=r'{0}'.format(legends[i_specie]))
    ax.legend(fontsize=6)

    fig.savefig(fig_path,
            dpi=300)
    
    plt.close()
    return None


def plot_min_max_avg(timestep, trace, fig_path, ylabel, title):
    """
    Description: Plot the min, average and max of single replciates among different replicates
    """
    fig_size = [87,87*l_w_ratio]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'time [s]',
                    fontsize=7,
                    labelpad=1.5)
    ax.set_ylabel(r''+ylabel.replace('_','_'),
                fontsize=7,
                labelpad=1.5)
    ax.set_title(r''+title.replace('_','_'), 
                fontsize=8,
                pad=4)
    ax.set_xlim(timestep[0] - 0.05*timestep[-1],
                timestep[-1])

    ax.set_ylim(0.9*np.amin(trace),
                1.2*np.amax(trace))
    
    tick_length = 2.0
    tick_width = 1.0
    ax.tick_params(labelsize=5,
                    length=tick_length,
                    width=tick_width,
                    direction='in',
                    left=True,
                    right=True,
                    bottom=True,
                    top=True,
                    which='major')


    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)

    cmap, c_space = choose_colormap(3)

    trace_min = np.min(trace, axis=1); trace_avg = np.mean(trace, axis=1); trace_max = np.max(trace, axis=1)
    
    min_avg_max= [trace_min, trace_avg, trace_max]
    
    legends = ['min', 'avg', 'max']
    for i in range(3):

        ax.plot(timestep,min_avg_max[i],
                            alpha=0.5,
                            linewidth=0.5,
                            zorder=2,
                            color=cmap(c_space[i]),
                            label=r'{0}'.format(legends[i]))
    
    ax.fill_between(timestep, trace_min, trace_max, 
                    color='lightblue', alpha=0.5,
                     zorder=1)
    
    ax.legend(fontsize=6)

    fig.savefig(fig_path,
            dpi=300)
    

    return None



def plot_histogram(data, figure_path, bins, xlabel, title):
    """
    Description: Plot the histogram of 1D numpy array data
    
    """
    fig_size = [87,87/1.618]

    fig = plt.figure(figsize=(fig_size[0]/25.4,fig_size[1]/25.4))

    plt.hist(data, bins=bins, alpha=0.7, color='limegreen')

    ax = plt.gca()

    xlabel = xlabel.replace('_','\_')
    ax.set_xlabel(r'{0}'.format(xlabel),
                  fontsize=7,
                  labelpad=1.5)

    ax.set_ylabel(r'{0}'.format('Frequency'),
                fontsize=7,
                labelpad=1.5)

    title = title.replace('_','\_')
    ax.set_title(r'{0}'.format(title),
                 fontsize=8,
                 pad=4)

    tick_length = 4.0
    tick_width = 1.5
    ax.tick_params(labelsize=5,
                    length=tick_length,
                    width=tick_width,
                    direction='in',
                    left=True,
                    right=True,
                    bottom=True,
                    top=True,
                    which='major')

    ax.tick_params(labelsize=5,
                    length=tick_length/1.5,
                    width=tick_width/1.5,
                    direction='in',
                    left=True,
                    right=False,
                    bottom=True,
                    top=False,
                    which='minor')

    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)

    mean = np.mean(data)
    plt.axvline(mean, color='red', linestyle='dashed', linewidth=1.5, label = 'Mean: {0:.3f}'.format(mean))

    median = np.median(data)
    plt.axvline(median, color='black', linestyle='dashed', linewidth=1.5, label = 'Median: {0:.3f}'.format(median))      

    min = np.min(data)
    plt.axvline(min, color='blue', linestyle='dashed', linewidth=1.5, label = 'Min: {0:.3f}'.format(min))  

    ax.legend(fontsize = 8)

    fig.savefig(figure_path, dpi = 300)

    return None