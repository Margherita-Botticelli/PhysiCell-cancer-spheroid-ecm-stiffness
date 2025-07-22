from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn

def plots_spheroid_area_growth_over_time(data, save_folder, title=True):
    """
    Plots the spheroid area growth over time, relative to the initial spheroid area
    
    Parameters:
    data (DataFrame): The input data containing simulation results.
    save_folder (str): Directory path where the plot will be saved.
    title (bool): Whether to include a title in the plot
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

    ecm_sensitivity = data['ecm_sensitivity'].iloc[0]
    chemotaxis_bias = data['chemotaxis_bias'].iloc[0]

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
    plt.figure(figsize=(6, 4), num=simulation)
    seaborn.set_context("paper")
    seaborn.set_style('ticks')
    plt.rcParams.update({'font.weight': 'bold',
        'axes.labelweight': 'bold'})

    #### Calculate mean and standard deviation of spheroid area ratio
    spheroid_area_mean = np.mean(spheroid_area_ratio, axis=0)
    # spheroid_area_std = np.std(spheroid_area_ratio, axis=0)
    spheroid_area_25_percentile = np.percentile(spheroid_area_ratio, 25 , axis=0) 
    spheroid_area_75_percentile = np.percentile(spheroid_area_ratio, 75 , axis=0)

    #### Plot mean and shaded area for standard deviation
    plt.plot(t, spheroid_area_mean, label=f'{ribose} mM', color=color_rib,lw=2)

    # spheroid_area_50_percentile = np.percentile(cell_count, 50 , axis=0)
    # plt.plot(t,spheroid_area_50_percentile,label=f'{ribose} mM',color='black',lw=2)
    
    # plt.fill_between(t, spheroid_area_mean + spheroid_area_std, spheroid_area_mean - spheroid_area_std, facecolor=color_rib, alpha=0.4)
    plt.fill_between(t, spheroid_area_75_percentile, spheroid_area_25_percentile,alpha=0.4, facecolor=color_rib)

    #### Set axis labels and ticks
    plt.ylabel(f'Growth relative to t$_0$', fontsize=15)
    plt.xlabel("Time [h]", fontsize=15)
    plt.yticks(np.arange(0, 8.1, 1), fontsize=15)
    plt.xticks(np.arange(0, t[-1] + 1, t[-1] / 4), fontsize=15)

    ax = plt.gca()  # Get current axis
    ax.spines["bottom"].set_linewidth(2)  # X-axis
    ax.spines["left"].set_linewidth(2)    # Y-axis
    ax.tick_params(axis='both', which='major', width=2)
    seaborn.despine()


    #### Set plot title
    if title:
        # plt.title(r'$\bf{Spheroid\,growth\,relative\,to\,t_0}$' + f'\n\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}', fontsize=12)
        plt.title(r'$\bf{Spheroid\,growth\,relative\,to\,t_0}$' + f'\nchemo={chemotaxis_bias}, ecm_sens={ecm_sensitivity}, s_mot={max_mot_speed}')

    
    plt.legend(fontsize=15, frameon=False, loc='upper left')

    #### Save the plot
    plt.savefig(save_folder + f'plots/spheroid_area_growth_over_time_{simulation}.png', bbox_inches="tight")


def plots_delaunay_mean_distance_over_time(data, save_folder, title=True):
    """
    Plots the mean Delaunay distance over time
    
    Parameters:
    data (DataFrame): The input data containing simulation results.
    save_folder (str): Directory path where the plot will be saved.
    title (bool): Whether to include a title in the plot
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

    
    ecm_sensitivity = data['ecm_sensitivity'].iloc[0]
    chemotaxis_bias = data['chemotaxis_bias'].iloc[0]


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
    plt.figure(figsize=(6,4), num=simulation + 1)
    seaborn.set_context("paper")
    seaborn.set_style('ticks')
    plt.rcParams.update({'font.weight': 'bold',
        'axes.labelweight': 'bold'})

    #### Calculate mean and standard deviation of Delaunay distance
    delaunay_distance_mean = np.mean(delaunay_distance, axis=0)
    # delaunay_distance_std = np.std(delaunay_distance, axis=0)
    delaunay_distance_25_percentile = np.percentile(delaunay_distance, 25 , axis=0) 
    delaunay_distance_75_percentile = np.percentile(delaunay_distance, 75 , axis=0)

    #### Plot mean and shaded area for standard deviation
    plt.plot(t, delaunay_distance_mean, label=f'{ribose} mM', color=color_rib,lw=2)

    # delaunay_distance_50_percentile = np.percentile(cell_count, 50 , axis=0)
    # plt.plot(t,delaunay_distance_50_percentile,label=f'{ribose} mM',color='black',lw=2)

    # plt.fill_between(t, delaunay_distance_mean + delaunay_distance_std, delaunay_distance_mean - delaunay_distance_std, facecolor=color_rib, alpha=0.4)
    plt.fill_between(t, delaunay_distance_75_percentile, delaunay_distance_25_percentile, facecolor=color_rib, alpha=0.4)

    #### Set axis labels and ticks
    plt.ylabel(r'Delaunay mean distance [$\mu$m]', fontsize=15)
    plt.xlabel("Time [h]", fontsize=15)
    plt.yticks(np.arange(15, 25.1, 5), fontsize=15)
    plt.xticks(np.arange(0, t[-1] + 1, t[-1] / 4), fontsize=15)

    ax = plt.gca()  # Get current axis
    ax.spines["bottom"].set_linewidth(2)  # X-axis
    ax.spines["left"].set_linewidth(2)    # Y-axis
    ax.tick_params(axis='both', which='major', width=2)
    seaborn.despine()


    #### Set plot title
    if title:
        # plt.title(r'$\bf{Delaunay\,mean\,distance}$' + f'\n\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}', fontsize=12)
        plt.title(r'$\bf{Delaunay\,mean\,distance}$' + f'\nchemo={chemotaxis_bias}, ecm_sens={ecm_sensitivity}, s_mot={max_mot_speed}')


    plt.legend(fontsize=15, frameon=False, loc='upper left')

    #### Save the plot
    plt.savefig(save_folder + f'plots/delaunay_mean_distance_over_time_{simulation}.png', bbox_inches="tight")


def plots_cell_count_over_time(data, save_folder, title=True):
    """
    Plots the number of cells over time
    
    Parameters:
    data (DataFrame): The input data containing simulation results.
    save_folder (str): Directory path where the plot will be saved.
    title (bool): Whether to include a title in the plot
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

    
    ecm_sensitivity = data['ecm_sensitivity'].iloc[0]
    chemotaxis_bias = data['chemotaxis_bias'].iloc[0]


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

    # cell_count_std = np.std(cell_count, axis=0)
    cell_count_25_percentile = np.percentile(cell_count, 25 , axis=0)
    cell_count_75_percentile = np.percentile(cell_count, 75 , axis=0)

    #### Create plot
    plt.figure(figsize=(6,4), num=simulation + 2)

    seaborn.set_context("paper")
    seaborn.set_style('ticks')
    plt.rcParams.update({'font.weight': 'bold',
        'axes.labelweight': 'bold'})

    #### Plot mean and shaded area for standard deviation
    plt.plot(t,cell_count_mean,label=f'{ribose} mM',color=color_rib,lw=2)

    # cell_count_50_percentile = np.percentile(cell_count, 50 , axis=0)
    # plt.plot(t,cell_count_50_percentile,label=f'{ribose} mM',color='black',lw=2)

    # plt.fill_between(t, cell_count_mean+cell_count_std, cell_count_mean-cell_count_std, facecolor=color_rib, alpha=0.4)
    plt.fill_between(t, cell_count_75_percentile, cell_count_25_percentile, alpha=0.4)
    
    #### Set axis labels and ticks
    plt.ylabel(f'Cell count', fontsize=15)
    plt.xlabel("Time [h]", fontsize=15)
    # plt.yticks(np.arange(0,1200,200))
    plt.yticks(fontsize=15)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/4), fontsize=15)

    ax = plt.gca()  # Get current axis
    ax.spines["bottom"].set_linewidth(2)  # X-axis
    ax.spines["left"].set_linewidth(2)    # Y-axis
    ax.tick_params(axis='both', which='major', width=2)
    seaborn.despine()


    #### Set plot title
    if title == True:
        # plt.title(r'$\bf{Cell\,number}$'+f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}', fontsize = 12)
        plt.title(r'$\bf{Cell\,number}$' + f'\nchemo={chemotaxis_bias}, ecm_sens={ecm_sensitivity}, s_mot={max_mot_speed}')

        
    plt.legend(fontsize=15, frameon=False, loc='upper left')
    
    #### Save the plot
    plt.savefig(save_folder + f'plots/cell_count_over_time_{simulation}.png', bbox_inches = "tight")



def plots_invasion_over_time(data, save_folder, title=True):
    """
    Plots of invasion over time
    
    Parameters:
    data (DataFrame): The input data containing simulation results.
    save_folder (str): Directory path where the plot will be saved.
    title (bool): Whether to include a title in the plot
    """
    #### Extract relevant parameters from the data
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = data['t'].unique() / 60  # Time in hours
    # prolif = data['prolif'].iloc[0]
    # cell_adh = data['cell_adh'].iloc[0]
    # cell_rep = data['cell_rep'].iloc[0]
    seeds = data['seed'].unique()
    # r_density = data['ecm_density_rate'].iloc[0]

    ecm_sensitivity = data['ecm_sensitivity'].iloc[0]
    chemotaxis_bias = data['chemotaxis_bias'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]


    #### Set color based on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    else: 
        color_rib = seaborn.color_palette('colorblind')[2]

    #### Collect cell count data
    invasion = []
    for seed in seeds:
        for timepoint in data['t'].unique(): 
            df_seed = data[(data['t'] == timepoint) & (data['seed'] == seed)]

            position_x = df_seed['position_x'].to_numpy()
            position_y = df_seed['position_y'].to_numpy()
            distances = np.sqrt(position_x**2 + position_y**2)
            max_distance = max(distances)
            invasion.append(max_distance)

    #### Reshape cell count data
    invasion = np.reshape(invasion, (-1, len(t))).astype(int)

    #### Calculate mean and standard deviation of cell count
    invasion_mean = np.mean(invasion, axis=0)
    invasion_25_percentile = np.percentile(invasion, 25 , axis=0)
    invasion_75_percentile = np.percentile(invasion, 75 , axis=0)

    #### Create plot
    plt.figure(figsize=(6,4), num=simulation + 3)

    seaborn.set_context("paper")
    seaborn.set_style('ticks')
    plt.rcParams.update({'font.weight': 'bold',
        'axes.labelweight': 'bold'})

    #### Plot mean and shaded area for standard deviation
    plt.plot(t,invasion_mean)
    # plt.fill_between(t, invasion_mean+invasion_std, invasion_mean-invasion_std, facecolor=color_rib, alpha=0.4)
    plt.fill_between(t, invasion_75_percentile, invasion_25_percentile, alpha=0.4)
    
    #### Set axis labels and ticks
    plt.ylabel(r'Invasion [$\mu$m]', fontsize=15)
    plt.xlabel("Time [h]", fontsize=15)
    plt.yticks(np.arange(0,801,200), fontsize=15)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/4), fontsize=15)

    ax = plt.gca()  # Get current axis
    ax.spines["bottom"].set_linewidth(2)  # X-axis
    ax.spines["left"].set_linewidth(2)    # Y-axis
    ax.tick_params(axis='both', which='major', width=2)
    seaborn.despine()


    #### Set plot title
    if title == True:
        plt.title(r'$\bf{Invasion}$' + f'\nchemo={chemotaxis_bias}, ecm_sens={ecm_sensitivity}, s_mot={max_mot_speed}')

        
    plt.legend(fontsize=15, frameon=False, loc='upper left')
    
    #### Save the plot
    plt.savefig(save_folder + f'plots/invasion_over_time_{simulation}.png', bbox_inches = "tight")
