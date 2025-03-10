import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import numpy as np
import matplotlib as mpl
import scipy.ndimage

def least_squares_function_approximated(data, save_folder):

    t = np.unique(data['t']).astype(int)

    #### Unique simulations in the dataset
    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))

    #### time values for experiments ribose 0 mM
    #### Invasive
    experiment = 'invasive'
    y_24 = 1.63
    y_48 = 2.65
    y_72 = 4.01
    y_96 = 6.00


    #### Non-invasive
    # experiment = 'noninvasive'
    # y_24 = 1.2
    # y_48 = 1.3
    # y_72 = 1.4
    # y_96 = 1.5

    max_mot_speed = []
    r_density = []
    least_squares = []

    for simulation in simulations:
        # print(simulation,flush=True)

        df_sim = data[data['simulation'] == simulation]
        # print(df_sim,flush=True)

        spheroid_area_init = df_sim[df_sim['t'] == min(t)]['spheroid_area'].to_numpy()

        spheroid_area_24 = df_sim[df_sim['t'] == 24*60]['spheroid_area'].to_numpy()
        spheroid_area_ratio_24 = spheroid_area_24 / spheroid_area_init
        f_24 = np.mean(spheroid_area_ratio_24)

        spheroid_area_48 = df_sim[df_sim['t'] == 48*60]['spheroid_area'].to_numpy()
        spheroid_area_ratio_48 = spheroid_area_48 / spheroid_area_init
        f_48 = np.mean(spheroid_area_ratio_48)

        spheroid_area_72 = df_sim[df_sim['t'] == 72*60]['spheroid_area'].to_numpy()
        spheroid_area_ratio_72 = spheroid_area_72 / spheroid_area_init
        f_72 = np.mean(spheroid_area_ratio_72)

        spheroid_area_96 = df_sim[df_sim['t'] == 96*60]['spheroid_area'].to_numpy()
        spheroid_area_ratio_96 = spheroid_area_96 / spheroid_area_init
        f_96 = np.mean(spheroid_area_ratio_96)

        #### least squares
        least_squares_sum = (y_24 - f_24)**2 + (y_48 - f_48)**2 + (y_72 - f_72)**2 + (y_96 - f_96)**2
        least_squares.append( least_squares_sum )

        print(f'{simulation=}',flush=True)
        print(y_24 ,f_24,flush=True)
        print(y_48 ,f_48,flush=True)
        print(y_72,f_72,flush=True)
        print(y_96 ,f_96,flush=True)
        print(least_squares_sum,flush=True)
        print('',flush=True)

        max_mot_speed.append(round(float(df_sim['max_mot_speed'].iloc[0]), 3))
        r_density.append(round(float(df_sim['ecm_density_rate'].iloc[0]), 5))

        # print(f'simulation={simulation}',flush=True)
        # print(f'least squares={least_squares}',flush=True)
        # print('\n',flush=True)

    fig, ax = plt.subplots(figsize=(7.5,5))
    seaborn.set_context("paper")
    seaborn.set_style('ticks')

    #### Create DataFrame for heatmap
    columns = np.unique(max_mot_speed)
    index = np.flip(np.unique(r_density))
    df = pd.DataFrame(columns=columns, index=index).fillna(0.0)
    annot_df = pd.DataFrame(columns=columns, index=index).fillna('NaN')

    for mot, degr, ls, sim in zip(max_mot_speed, r_density, least_squares,simulations):
        df[mot][degr] = ls
        annot = f"{sim}"
        annot_df[mot][degr] = annot

    annot_arr = annot_df.to_numpy()


    color_rib = seaborn.color_palette('colorblind')[0]
    color_rib_dark = seaborn.color_palette('dark')[0]

    #### Define colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white', color_rib, color_rib_dark])

    #### Plot heatmap
    hmap = seaborn.heatmap(df,cmap=cmap,ax=ax,vmax=0.2, annot=annot_arr, fmt="s")

    hmap.figure.axes[-1].yaxis.set_label_position('left')
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(label='Least squares',fontsize=12)

    
    #### Set axis labels and title
    plt.ylabel(r'$r_{deg,0}$ [min$^{-1}$]', color='black', fontsize=12)
    plt.yticks(color='black', fontsize=12,rotation=45,va='top')
    plt.xlabel(r'$S_0$ [$\mu$m$\cdot$min$^{-1}$]', color='black', fontsize=12)
    plt.xticks(color='black', fontsize=12)

    plt.savefig(save_folder + f'plots/least_squares_heatmap_{experiment}.png', bbox_inches="tight")
    plt.close()



def least_squares_function_cell_count(data, save_folder):

    t = np.unique(data['t']).astype(int)

    #### Unique simulations in the dataset
    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))

    #### time values for experiments ribose 0 mM
    #### Invasive
    experiment = 'invasive_cell_count'

    cell_values_24 = [199,246,280,209,188,170,170,219,212,226,239,236]
    cell_values_24 = np.array(cell_values_24)
    y_24 = np.mean(cell_values_24)

    cell_values_72 = [462,543,487,601,490,619,365,614,321,269,300,490,271,357,214,305,769,1167,950,1098,725,967,699,863]
    cell_values_72 = np.array(cell_values_72)
    y_72 = np.mean(cell_values_72)

    #### Non-invasive
    # experiment = 'noninvasive_cell_count'

    max_mot_speed = []
    r_density = []
    least_squares = []

    for simulation in simulations:
        # print(simulation,flush=True)
        # print(data[(data['simulation'] == simulation)],flush=True)
        df_sim = data[data['simulation'] == simulation]

        df_sim_24 = data[(data['simulation'] == simulation) & (data['t'] == 24*60)]
        cell_count_24 = df_sim_24['ID'].iloc[-1] +1
        f_24 = np.mean(cell_count_24)

        df_sim_72 = data[(data['simulation'] == simulation) & (data['t'] == 72*60)]
        cell_count_72 = df_sim_72['ID'].iloc[-1] + 1
        f_72 = np.mean(cell_count_72)

        #### least squares
        least_squares_sum = (y_24 - f_24)**2 + (y_72 - f_72)**2 
        least_squares.append( least_squares_sum )

        # print(f'{simulation=}',flush=True)
        # print(y_24 ,f_24,flush=True)
        # print(y_72,f_72,flush=True)
        # print(least_squares_sum,flush=True)
        # print('',flush=True)

        max_mot_speed.append(round(float(df_sim['max_mot_speed'].iloc[0]), 3))
        r_density.append(round(float(df_sim['ecm_density_rate'].iloc[0]), 5))

        # print(f'simulation={simulation}',flush=True)
        # print(f'least squares={least_squares}',flush=True)
        # print('\n',flush=True)

    fig, ax = plt.subplots(figsize=(7.5,5))
    seaborn.set_context("paper")
    seaborn.set_style('ticks')

    #### Create DataFrame for heatmap
    columns = np.unique(max_mot_speed)
    index = np.flip(np.unique(r_density))
    df = pd.DataFrame(columns=columns, index=index).fillna(0.0)
    annot_df = pd.DataFrame(columns=columns, index=index).fillna('NaN')

    for mot, degr, ls, sim in zip(max_mot_speed, r_density, least_squares,simulations):
        df[mot][degr] = ls
        annot = f"{sim}"
        annot_df[mot][degr] = annot

    annot_arr = annot_df.to_numpy()


    color_rib = seaborn.color_palette('colorblind')[0]
    color_rib_dark = seaborn.color_palette('dark')[0]

    #### Define colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white', color_rib, color_rib_dark])

    #### Plot heatmap
    hmap = seaborn.heatmap(df,cmap=cmap,ax=ax,annot=annot_arr, fmt="s",vmax=10000)

    hmap.figure.axes[-1].yaxis.set_label_position('left')
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(label='Least squares',fontsize=12)

    
    #### Set axis labels and title
    plt.ylabel(r'$r_{deg,0}$ [min$^{-1}$]', color='black', fontsize=12)
    plt.yticks(color='black', fontsize=12,rotation=45,va='top')
    plt.xlabel(r'$S_0$ [$\mu$m$\cdot$min$^{-1}$]', color='black', fontsize=12)
    plt.xticks(color='black', fontsize=12)

    plt.savefig(save_folder + f'plots/least_squares_heatmap_{experiment}.png', bbox_inches="tight")
    plt.close()