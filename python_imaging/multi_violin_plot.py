from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import splrep, splev
import seaborn as sns
import pandas as pd
import sys

def cell_number_data(snapshot, data_folder):
    # Function to get data for the plot
    
    #################################  Load data  ########################
    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', data_folder)

    # loads and reads ECM data
    mcds.load_ecm(snapshot + '_ECM.mat', data_folder)

    # time
    t = ( mcds.get_time() )

    # Cells
    cell_df = mcds.get_cell_df()
    
    # Get each cell's position and total speed
    positions_x = []
    positions_y = []
    total_speeds = []
    dist = []
    
    for j in cell_df.index:
        # Get data from cell
        position_x = cell_df.loc[j, 'position_x']
        position_y = cell_df.loc[j, 'position_y']
        dist_temp = np.sqrt(position_x**2 + position_y**2)
        
        # Append data to lists
        dist.append(dist_temp)
       
    return t, dist

def cell_number_plot(snapshot,ribose,simulation):
    # Violin plot of distance from origin of each cell for different ribose concentrations
    
    # Initiate distance array
    n = len(ribose)
    dist =  [[] for _ in range(n)]
    
    for i in np.arange(0,n):
        
        # Data folder
        data_folder = '../ribose_' + ribose[i] + '/output/output_' + simulation  
        
        # Save folder
        save_folder = '../results/'    
        
        # Call data function
        t, dist[i] = cell_number_data(snapshot, data_folder)
        
        
    # distance =  [item for sublist in dist for item in sublist]
    
    # Build data disctionary
    data = {'ribose':  [],
            'distance': []}
    
    # Append data for plot
    for r, dis in zip(ribose, dist):
        for d in dis:
            data['ribose'].append(r)
            data['distance'].append(d)
      
    # Build data frame
    df = pd.DataFrame(data)
    
    # Initiate figure
    fig, ax = plt.subplots(figsize=(9, 6))
    
    # Produce violin plot
    sns.violinplot(x='ribose', y='distance', data=df, cut=0,
                   scale='width', inner=None, linewidth=1, 
                   saturation=1)
    
    # Add dots for each cell on top
    sns.swarmplot(x='ribose', y='distance', data=df, linewidth=1,color='black')
    
    # Ticks
    _ = plt.xticks(fontsize=15)
    sns.despine(left=True)
    plt.yticks(np.arange(0,151,20),fontsize=15)

    # Axis labels
    ax.set_ylabel('Distance [micron]',fontsize=20)
    ax.set_xlabel('Ribose [mM]',fontsize=20)
    
    # Plot title
    plt.title('Distribution of cells from origin at t={t}min'.format(t=t),fontsize=20)
        
    plt.savefig(save_folder + 'multi_violin_plot_' + simulation + '.png', bbox_inches = "tight")
    

if __name__ == '__main__':
    
    # Figure resolution
    mpl.rcParams['figure.dpi'] = 300
    
    # Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')
    
    # Snapshot 
    snapshot = 'output00000300'
    
    # Ribose array
    # ribose = ['0','25','50','100','150','200']
    # ribose = ['0','50']
    ribose = ['0','50','200']
    
    # Simulation
    simulation = '9_0'
    
    # Distance violin plot
    cell_number_plot(snapshot,ribose,simulation)
    
    # Make background transparent
    # fig = plt.gcf()
    # fig.patch.set_alpha(0.0)

        