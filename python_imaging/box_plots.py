from cgitb import small
from pyMCDS_ECM import *
import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import numpy as np
import matplotlib as mpl

def box_plot_spheroid_area_growth(data, save_folder, title=True):
    """
    Boxplot of spheroid area growth relative to initial time.
    
    Parameters:
    data (DataFrame): The input data containing simulation results.
    save_folder (str): Directory path where the plot will be saved.
    title (bool): Whether to include a title in the plot.
    """
    #### Filter data to include only rows where cell ID is 0
    data = data[(data['ID'] == 0)]

    #### Extract relevant parameters from the data
    t = data['t'].unique() / 60  # Time in hours
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    r_density = data['ecm_density_rate'].iloc[0]

    #### Get values of ribose concentrations and simulations
    riboses = list(dict.fromkeys(data['ribose'].values.tolist()))
    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))

    #### Prepare a list to store reshaped data
    df_list = []
    
    for simulation in simulations:
        for ribose in riboses:
            #### Filter data for each ribose concentration and simulation
            df = data[(data['ribose'] == ribose) & (data['simulation'] == simulation)]
            spheroid_area = df['spheroid_area'].to_numpy()

            #### Reshape spheroid area data into a 2D array with timepoints as columns
            spheroid_area = np.reshape(spheroid_area, (-1, len(t))).astype(float)
            spheroid_area_ratio = spheroid_area[:, 0]
            spheroid_area_ratio = np.reshape(spheroid_area_ratio, (-1, 1))
            spheroid_area_ratio = spheroid_area / spheroid_area_ratio
            spheroid_area_ratio = np.array(spheroid_area_ratio)

            df1 = pd.DataFrame() # Temporary DataFrame to store ratios

            #### Add each time point as a separate column
            for i in range(spheroid_area_ratio.shape[1]):
                df1[f'{i}'] = spheroid_area_ratio[:,i]

            #### Include simulation and ribose values as columns for identification
            df1['simulation'] = simulation
            df1['ribose'] = ribose

            #### Append the processed DataFrame to the list
            df_list.append(df1)

    #### Concatenate all DataFrames into a single DataFrame
    df = pd.concat(df_list, copy=False, axis=0)

    #### Create plot
    fig = plt.figure(figsize=(4, 4), num=simulation)
    plt.rcParams.update({'font.weight': 'bold',
        'axes.labelweight': 'bold'})

    #### Set seaborn context and style
    seaborn.set_context("paper")
    seaborn.set_style('ticks')
    # seaborn.despine(top=True,right=True)
    seaborn.despine()


    #### Melt DataFrame
    df = pd.melt(df, id_vars=['simulation', 'ribose'], var_name='t', value_name='spheroid_area_ratio', ignore_index=True)

    #### Filter for a specific time point
    df = df[df['t'].isin(['72'])]

    #### Create a boxplot of the spheroid area growth ratios
    meanlineprops = dict(linestyle='-',linewidth=4)
    box_plot = seaborn.boxplot(data=df, x='ribose', y='spheroid_area_ratio', hue='simulation', fliersize=0,legend=False, linewidth=4, gap=0.2, showmeans=True, meanline=True, meanprops=meanlineprops)

    #### Define colors and hatches for the box plots
    color_rib_0 = seaborn.color_palette('colorblind')[0]
    color_rib_50 = seaborn.color_palette('colorblind')[1]
    color_rib_200 = seaborn.color_palette('colorblind')[2]

    # color_t_24 = seaborn.color_palette('colorblind')[3]
    # color_t_72 = seaborn.color_palette('colorblind')[4]

    # palette = [color_t_72,color_t_72, color_t_72, color_t_72,color_t_72, color_t_72]

    palette = [color_rib_0, color_rib_50, color_rib_200,color_rib_0,color_rib_50,color_rib_200]
    # hatches = [ '\\\\\\\\\\', '\\\\\\\\\\','\\\\\\\\\\', '\\\\\\\\\\','\\\\\\\\\\', '\\\\\\\\\\']

    #### Customize box plot appearance based on colors and hatches
    for i,box in enumerate(box_plot.patches):
        color = palette[i]
        # hatch = hatches[i]
        box.set_edgecolor(color)
        box.set_facecolor(mpl.colors.to_rgba(color, 0.5)) # Semi-transparent fill
        # box.set_hatch(hatch)

        #### Iterate over whiskers and median lines
        for j in range(7*i,7*(i+1)):
            if j%(7*i+4)==0 and j!= 0:
                box_plot.lines[j].set_color(mpl.colors.to_rgba(color,0))
            else:
                box_plot.lines[j].set_color(color)

    #### Create a custom legend for ribose concentrations
    legend_elements = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_rib_0,0.5), edgecolor=color_rib_0,  lw=4, label='0 mM'), 
            mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_rib_50,0.5),edgecolor=color_rib_50, lw=4, label='50 mM'), 
            mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_rib_200,0.5),edgecolor=color_rib_200, lw=4, label='200 mM')]

    # legend_elements = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_t_72,0.5),edgecolor=color_t_72, lw=4, label='72 h')]

    plt.legend(handles=legend_elements,frameon=False, loc='upper right',fontsize=15)

    #### Set axis labels and ticks
    plt.ylabel(f'Growth relative to t$_0$',fontsize=15)
    plt.yticks(np.arange(0, 6.1, 2),fontsize=15)
    box_plot.set_xticks([-0.2,0.2,0.8,1.2,1.8,2.2]) 
    box_plot.set_xticklabels(['Control', 'GM6001','Control', 'GM6001','Control', 'GM6001'], rotation=45,ha='right',fontsize=15) 

    ax = plt.gca()  # Get current axis
    ax.spines["bottom"].set_linewidth(2)  # X-axis
    ax.spines["left"].set_linewidth(2)    # Y-axis
    ax.tick_params(axis='both', which='major', width=2)
    
    seaborn.despine()

    #### Set plot title
    if title:
        plt.title(r'$\bf{Spheroid\,growth\,relative\,to\,t_0}$' + f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}')
    
    #### Save the plot
    plt.savefig(save_folder + f'plots/spheroid_area_growth_box_plot_{simulation}.png', bbox_inches="tight",dpi=600)



def box_plot_delaunay_mean_distance(data, save_folder, title=True):
    """
    Boxplot of Delaunay mean distance
    
    Parameters:
    data (DataFrame): The input data containing simulation results.
    save_folder (str): Directory path where the plot will be saved.
    title (bool): Whether to include a title in the plot
    """
    #### Filter data for cell ID == 0
    data = data[(data['ID'] == 0)]

    #### Extract relevant parameters from the data
    t = data['t'].unique() / 60  # Time in hours
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    r_density = data['ecm_density_rate'].iloc[0]

    #### Get values of ribose concentrations and simulations
    riboses = list(dict.fromkeys(data['ribose'].values.tolist()))
    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))

    #### Prepare a list to store reshaped data
    df_list = []
    
    for simulation in simulations:
        for ribose in riboses:
          
            #### Filter data for each ribose concentration and simulation
            df = data[(data['ribose'] == ribose) & (data['simulation'] == simulation)]
            delaunay_distance = df['delaunay_distance'].to_numpy()

            #### Reshape delaunay distance data into a 2D array with timepoints as columns
            delaunay_distance = np.reshape(delaunay_distance, (-1, len(t))).astype(float)
            delaunay_distance = np.array(delaunay_distance)

            df1 = pd.DataFrame() # Temporary DataFrame

            #### Add each time point as a separate column
            for i in range(delaunay_distance.shape[1]):
                df1[f'{i}'] = delaunay_distance[:,i]

            #### Include simulation and ribose values as columns for identification
            df1['simulation'] = simulation
            df1['ribose'] = ribose

            #### Append the processed DataFrame to the list
            df_list.append(df1)

    #### Concatenate all DataFrames into a single DataFrame
    df = pd.concat(df_list, copy=False, axis=0)

    #### Create plot
    fig = plt.figure(figsize=(4, 4), num=simulation+1)
    plt.rcParams.update({'font.weight': 'bold',
        'axes.labelweight': 'bold'})

    #### Set seaborn context and style
    seaborn.set_context("paper")
    seaborn.set_style('ticks')
    # seaborn.despine(top=True,right=True)
    seaborn.despine()

    #### Melt DataFrame
    df = pd.melt(df, id_vars=['simulation', 'ribose'], var_name='t', value_name='delaunay_distance', ignore_index=True)

    #### Filter for specific time points
    if simulation == 293:
        df = df[df['t'].isin(['24','72'])]
    elif simulation == 240:
        df = df[df['t'].isin(['24','96'])]

    #### Create a boxplot of the Delaunay mean distance
    meanlineprops = dict(linestyle='-',linewidth=4)
    box_plot = seaborn.boxplot(data=df, x='ribose', y='delaunay_distance', hue='t', fliersize=0,legend=False, linewidth=4,gap=0.2, showmeans=True, meanline=True, meanprops=meanlineprops)

    #### Define colors and hatches for the box plots
    color_rib_0 = seaborn.color_palette('colorblind')[0]
    color_rib_50 = seaborn.color_palette('colorblind')[1]
    color_rib_200 = seaborn.color_palette('colorblind')[2]

    color_t_24 = seaborn.color_palette('colorblind')[5]
    color_t_24 = seaborn.color_palette('colorblind')[3]
    color_t_72 = seaborn.color_palette('colorblind')[4]
    color_t_96 = seaborn.color_palette('colorblind')[9]

    if simulation == 293:
        palette = [color_t_24, color_t_24,color_t_24, color_t_72,color_t_72, color_t_72]
    elif simulation == 240:
        palette = [color_t_24, color_t_24,color_t_24, color_t_96,color_t_96, color_t_96]
    # hatches = ['','','','\\\\\\\\', '\\\\\\\\', '\\\\\\\\']
    # hatches = ['','','','//////', '//////', '//////']

    #### Customize box plot appearance based on colors and hatches
    for i,box in enumerate(box_plot.patches):
        color = palette[i]
        # hatch = hatches[i]
        box.set_edgecolor(color)
        box.set_facecolor(mpl.colors.to_rgba(color, 0.5))
        # box.set_hatch(hatch)
        #### Iterate over whiskers and median lines
        for j in range(7*i,7*(i+1)):
            if j%(7*i+4)==0 and j!= 0:
                box_plot.lines[j].set_color(mpl.colors.to_rgba(color,0))
            else:
                box_plot.lines[j].set_color(color)

    #### Create a custom legend for ribose concentrations
    # legend_elements = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_rib_0,0.5), edgecolor=color_rib_0,  lw=4, label='0 mM'), 
    #         mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_rib_50,0.5),edgecolor=color_rib_50, lw=4, label='50 mM'), 
    #         mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_rib_200,0.5),edgecolor=color_rib_200, lw=4, label='200 mM')]

    if simulation == 293:
        legend_elements = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_t_24,0.5), edgecolor=color_t_24,  lw=4, label='24 h'), 
            mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_t_72,0.5),edgecolor=color_t_72, lw=4, label='72 h')]
    elif simulation == 240:
        legend_elements = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_t_24,0.5), edgecolor=color_t_24,  lw=4, label='24 h'), 
            mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_t_96,0.5),edgecolor=color_t_96, lw=4, label='96 h')]

    plt.legend(handles=legend_elements, ncol=2, frameon=False, loc='upper center',bbox_to_anchor=(0.5, 1.2),fontsize=15)

    #### Set axis labels and ticks
    plt.ylabel(r'Delaunay mean distance [$\mu$m]',fontsize=15)
    if simulation == 240:
        plt.yticks(np.arange(12, 28.1, 4),fontsize=15)
    elif simulation == 293:
        plt.yticks(np.arange(15, 30.1, 5),fontsize=15)

    plt.xlabel("Ribose [mM]",fontsize=15)
    # box_plot.set_xticks([-0.2,0.2,0.8,1.2,1.8,2.2]) 
    box_plot.set_xticks([0,1,2]) 
    box_plot.set_xticklabels([0,50,200],fontsize=15) 
    # box_plot.set_xticklabels([24,72,24,72,24,72]) 
    # box_plot.set_xticklabels([24,96,24,96,24,96]) 

    ax = plt.gca()  # Get current axis
    ax.spines["bottom"].set_linewidth(2)  # X-axis
    ax.spines["left"].set_linewidth(2)    # Y-axis
    ax.tick_params(axis='both', which='major', width=2)
    
    seaborn.despine()

    #### Set plot title
    if title:
        plt.title(r'$\bf{Delaunay\,mean\,distance}$' + f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}')
        
    #### Save the plot
    plt.savefig(save_folder + f'plots/delaunay_mean_distance_box_plot_{simulation}.png', bbox_inches="tight")


def box_plot_cell_count(data, save_folder, title=True):
    """
    Boxplot of cell count

    Parameters:
    data (DataFrame): The input data containing simulation results.
    save_folder (str): Directory path where the plot will be saved.
    title (bool): Whether to include a title in the plot
    """

    #### Extract relevant parameters from the data
    t = data['t'].unique() / 60  # Time in hours
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    r_density = data['ecm_density_rate'].iloc[0]
       
    #### Get values of ribose concentrations and simulations
    riboses = list(dict.fromkeys(data['ribose'].values.tolist()))
    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))

    #### Prepare a list to store reshaped data
    df_list = []
    
    for simulation in simulations:
        for ribose in riboses:
          
            #### Filter data for each ribose concentration and simulation
            df = data[(data['ribose'] == ribose) & (data['simulation'] == simulation)]

            #### Get values of seeds
            seeds = df['seed'].unique()

            #### Collect cell count data
            cell_count = []
            for seed in seeds:
                for timepoint in df['t'].unique(): 
                    df_seed = df[(df['t'] == timepoint) & (df['seed'] == seed)]
                    cell_count_t = df_seed['ID'].iloc[-1]
                    cell_count.append(cell_count_t + 1)

            #### Reshape cell count data
            cell_count = np.reshape(cell_count, (-1, len(t))).astype(int)
            cell_count = np.array(cell_count)

            df1 = pd.DataFrame() # Temporary DataFrame

            #### Add each time point as a separate column
            for i in range(cell_count.shape[1]):
                df1[f'{i}'] = cell_count[:,i]

            #### Include simulation and ribose values as columns for identification
            df1['simulation'] = simulation
            df1['ribose'] = ribose

            #### Append the processed DataFrame to the list
            df_list.append(df1)

    #### Concatenate all DataFrames into a single DataFrame
    df = pd.concat(df_list, copy=False, axis=0)

    #### Create plot
    fig = plt.figure(figsize=(4, 4), num=simulation+2)
    plt.rcParams.update({'font.weight': 'bold',
        'axes.labelweight': 'bold'})

    #### Set seaborn context and style
    seaborn.set_context("paper")
    seaborn.set_style('ticks')

    #### Melt DataFrame
    df = pd.melt(df, id_vars=['simulation', 'ribose'], var_name='t', value_name='cell_count', ignore_index=True)

    #### Filter for specific time points
    if simulation == 293:
        df = df[df['t'].isin(['24','72'])]
    elif simulation == 240:
        df = df[df['t'].isin(['24','96'])]

    #### Create a boxplot of the cell count
    meanlineprops = dict(linestyle='-',linewidth=4)
    box_plot = seaborn.boxplot(data=df, x='ribose', y='cell_count', hue='t', fliersize=0,legend=False, linewidth=4, gap=0.2, showmeans=True, meanline=True, meanprops=meanlineprops)

    #### Define colors and hatches for the box plots
    color_rib_0 = seaborn.color_palette('colorblind')[0]
    color_rib_50 = seaborn.color_palette('colorblind')[1]
    color_rib_200 = seaborn.color_palette('colorblind')[2]

    color_t_24 = seaborn.color_palette('colorblind')[3]
    color_t_72 = seaborn.color_palette('colorblind')[4]
    color_t_96 = seaborn.color_palette('colorblind')[9]

    if simulation == 293:
        palette = [color_t_24, color_t_24, color_t_24, color_t_72, color_t_72, color_t_72]
    elif simulation == 240:
        palette = [color_t_24, color_t_24, color_t_24, color_t_96, color_t_96, color_t_96]

    # palette = [color_rib_0, color_rib_50, color_rib_200,color_rib_0,color_rib_50,color_rib_200]
    # hatches = ['','','','\\\\\\\\', '\\\\\\\\', '\\\\\\\\']
    # hatches = ['','','','//////', '//////', '//////']

    #### Customize box plot appearance based on colors and hatches
    for i,box in enumerate(box_plot.patches):
        print(f'{i=}')
        color = palette[i]
        # hatch = hatches[i]
        box.set_edgecolor(color)
        box.set_facecolor(mpl.colors.to_rgba(color, 0.5))
        # box.set_hatch(hatch)

        #### Iterate over whiskers and median lines
        for j in range(7*i,7*(i+1)):
            if j%(7*i+4)==0 and j!= 0:
                box_plot.lines[j].set_color(mpl.colors.to_rgba(color,0))
            else:
                box_plot.lines[j].set_color(color)


    #### Create a custom legend for ribose concentrations
    # legend_elements = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_rib_0,0.5), edgecolor=color_rib_0,  lw=4, label='0 mM'), 
    #         mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_rib_50,0.5),edgecolor=color_rib_50, lw=4, label='50 mM'), 
    #         mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_rib_200,0.5),edgecolor=color_rib_200, lw=4, label='200 mM')]

    if simulation == 293:
        legend_elements = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_t_24,0.5), edgecolor=color_t_24,  lw=4, label='24 h'), 
            mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_t_72,0.5),edgecolor=color_t_72, lw=4, label='72 h')]
    elif simulation == 240:
        legend_elements = [mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_t_24,0.5), edgecolor=color_t_24,  lw=4, label='24 h'), 
            mpl.patches.Patch(facecolor=mpl.colors.to_rgba(color_t_96,0.5),edgecolor=color_t_96, lw=4, label='96 h')]

    plt.legend(handles=legend_elements, ncol=2, frameon=False, loc='upper center',bbox_to_anchor=(0.5, 1.2),fontsize=15)

    #### Set axis labels and ticks
    plt.ylabel(f'Cell count',fontsize=15)
    if simulation == 240:
        box_plot.set_yticks([0,240,300,500,700]) 
        box_plot.set_yticklabels(['',240,300,500,700],fontsize=15) 
    elif simulation == 293:
        plt.yticks(np.arange(0, 1500.1, 500),fontsize=15)
    plt.xlabel("Ribose [mM]",fontsize=15)
    box_plot.set_xticks([0,1,2]) 
    box_plot.set_xticklabels([0,50,200],fontsize=15) 
    # plt.xlabel("Time [h]")
    # box_plot.set_xticks([-0.2,0.2,0.8,1.2,1.8,2.2]) 
    # box_plot.set_xticklabels([24,72,24,72,24,72]) 
    # box_plot.set_xticklabels([24,96,24,96,24,96]) 

    ax = plt.gca()  # Get current axis
    ax.spines["bottom"].set_linewidth(2)  # X-axis
    ax.spines["left"].set_linewidth(2)    # Y-axis
    ax.tick_params(axis='both', which='major', width=2)

    seaborn.despine()

    #### Set plot title
    if title == True:
        plt.title(r'$\bf{Cell\,number}$'+f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}')

    #### Save the plot
    plt.savefig(save_folder + f'plots/cell_count_box_plot_{simulation}.png', bbox_inches = "tight")
