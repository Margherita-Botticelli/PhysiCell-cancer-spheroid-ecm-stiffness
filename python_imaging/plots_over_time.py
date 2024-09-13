from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn

def plots_spheroid_area_growth_over_time(data, save_folder, title=True):
    """
    Plots the spheroid area growth over time, relative to the initial spheroid area.
    
    Parameters:
    data (DataFrame): The input data containing simulation results.
    save_folder (str): Directory path where the plot will be saved.
    title (bool): Whether to include a title in the plot.
    """
    #### Filter data for cell ID == 0 
    data = data[(data['ID'] == 0)]

    #### Extract relevant parameters from the data
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = data['t'].unique() / 60  # Time in hours
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    spheroid_area = data['spheroid_area'].to_numpy()
    r_density = data['ecm_density_rate'].iloc[0]

    #### Set color based on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    else: 
        color_rib = seaborn.color_palette('colorblind')[2]

    #### Reshape spheroid area data
    spheroid_area = np.reshape(spheroid_area, (-1, len(t))).astype(float)
    spheroid_area_ratio = spheroid_area[:, 0]
    spheroid_area_ratio = np.reshape(spheroid_area_ratio, (-1, 1))
    spheroid_area_ratio = spheroid_area / spheroid_area_ratio

    #### Create plot
    plt.figure(figsize=(5, 3), num=simulation)

    #### Calculate mean and standard deviation of spheroid area ratio
    spheroid_area_mean = np.mean(spheroid_area_ratio, axis=0)
    spheroid_area_std = np.std(spheroid_area_ratio, axis=0)

    #### Plot mean and shaded area for standard deviation
    plt.plot(t, spheroid_area_mean, label=f'{ribose} mM', color=color_rib)
    plt.fill_between(t, spheroid_area_mean + spheroid_area_std, spheroid_area_mean - spheroid_area_std, facecolor=color_rib, alpha=0.4)

    #### Set axis labels and ticks
    plt.ylabel(f'Growth relative to t$_0$', fontsize=18)
    plt.xlabel("Time [h]", fontsize=18)
    plt.yticks(np.arange(0, 8.1, 1), fontsize=15)
    plt.xticks(np.arange(0, t[-1] + 1, t[-1] / 4), fontsize=15)

    #### Set plot title
    if title:
        plt.title(r'$\bf{Spheroid\,growth\,relative\,to\,t_0}$' + f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}', fontsize=12)
    
    plt.legend(title=f'Ribose', fontsize=15, title_fontsize=15)

    #### Save the plot
    plt.savefig(save_folder + f'plots/spheroid_area_growth_over_time_{simulation}.png', bbox_inches="tight")


def plots_delaunay_mean_distance_over_time(data, save_folder, title=True):
    """
    Plots the mean Delaunay distance over time.
    
    Parameters:
    data (DataFrame): The input data containing simulation results.
    save_folder (str): Directory path where the plot will be saved.
    title (bool): Whether to include a title in the plot.
    """
    #### Filter data for cell ID == 0
    data = data[(data['ID'] == 0)]
    
    #### Extract relevant parameters from the data
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = data['t'].unique() / 60  # Time in hours
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    delaunay_distance = data['delaunay_distance'].to_numpy()
    r_density = data['ecm_density_rate'].iloc[0]

    #### Set color based on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    else: 
        color_rib = seaborn.color_palette('colorblind')[2]

    #### Reshape Delaunay distance data
    delaunay_distance = np.reshape(delaunay_distance, (-1, len(t))).astype(float)

    #### Create plot
    plt.figure(figsize=(5, 4), num=simulation + 1)

    #### Calculate mean and standard deviation of Delaunay distance
    delaunay_distance_mean = np.mean(delaunay_distance, axis=0)
    delaunay_distance_std = np.std(delaunay_distance, axis=0)

    #### Plot mean and shaded area for standard deviation
    plt.plot(t, delaunay_distance_mean, label=f'{ribose} mM', color=color_rib)
    plt.fill_between(t, delaunay_distance_mean + delaunay_distance_std, delaunay_distance_mean - delaunay_distance_std, facecolor=color_rib, alpha=0.4)

    #### Set axis labels and ticks
    plt.ylabel(r'Delaunay mean distance [$\mu$m]', fontsize=18)
    plt.xlabel("Time [h]", fontsize=18)
    plt.yticks(np.arange(0, 25.1, 5), fontsize=15)
    plt.xticks(np.arange(0, t[-1] + 1, t[-1] / 4), fontsize=15)

    #### Set plot title
    if title:
        plt.title(r'$\bf{Delaunay\,mean\,distance}$' + f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}', fontsize=12)
        
    plt.legend(title=f'Ribose', fontsize=15, title_fontsize=15)

    #### Save the plot
    plt.savefig(save_folder + f'plots/delaunay_mean_distance_over_time_{simulation}.png', bbox_inches="tight")


def plots_cell_count_over_time(data, save_folder, title=True):
    """
    Plots the number of cells over time.
    
    Parameters:
    data (DataFrame): The input data containing simulation results.
    save_folder (str): Directory path where the plot will be saved.
    title (bool): Whether to include a title in the plot.
    """
    #### Extract relevant parameters from the data
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = data['t'].unique() / 60  # Time in hours
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    seeds = data['seed'].unique()
    r_density = data['ecm_density_rate'].iloc[0]

    #### Set color based on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    else: 
        color_rib = seaborn.color_palette('colorblind')[2]

    #### Collect cell count data
    cell_count = []
    for seed in seeds:
        for timepoint in data['t'].unique(): 
            df_seed = data[(data['t'] == timepoint) & (data['seed'] == seed)]
            cell_count_t = df_seed['ID'].iloc[-1]
            cell_count.append(cell_count_t + 1)

    #### Reshape cell count data
    cell_count = np.reshape(cell_count, (-1, len(t))).astype(int)

    #### Calculate mean and standard deviation of cell count
    cell_count_mean = np.mean(cell_count, axis=0)
    cell_count_std = np.std(cell_count, axis=0)

    #### Create plot
    plt.figure(figsize=(5, 4), num=simulation)

    #### Plot mean and shaded area for standard deviation
    plt.plot(t,cell_count_mean,label=f'{ribose} mM',color=color_rib)
    plt.fill_between(t, cell_count_mean+cell_count_std, cell_count_mean-cell_count_std, facecolor=color_rib, alpha=0.4)
    
    #### Set axis labels and ticks
    plt.ylabel(f'Cell count',fontsize=18)
    plt.xlabel("Time [h]",fontsize=18)
    plt.yticks(np.arange(0,1200,200),fontsize=15)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/4),fontsize=15)

    #### Set plot title
    if title == True:
        plt.title(r'$\bf{Cell\,number}$'+f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}', fontsize = 12)
        
    plt.legend(title=f'Ribose',fontsize=15,title_fontsize=15)
    
    #### Save the plot
    plt.savefig(save_folder + f'plots/cell_count_over_time_{simulation}.png', bbox_inches = "tight")

