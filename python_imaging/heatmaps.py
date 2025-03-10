import matplotlib.pyplot as plt
import seaborn
import pandas as pd
import numpy as np
import matplotlib as mpl
import scipy.ndimage

def plots_speed_vs_degr_delaunay(data, simulation_name, save_folder, title=True, contour=False):
    """
    Plot the relationship between max migration speed and degradation rate for the Delaunay distance

    Parameters:
    - data: pandas DataFrame containing the simulation data.
    - simulation_name: String identifier for the simulation.
    - save_folder: Directory path to save the plot.
    - title: Boolean to decide whether to include a title on the plot
    """
    # plt.figure()
    fig, ax = plt.subplots(figsize=(7.5,5))
    seaborn.set_context("paper")
    seaborn.set_style('ticks')

    #### Collect simulation parameters from the DataFrame
    ribose = data['ribose'].iloc[0]
    t = np.unique(data['t']).astype(int)
    prolif = round(float(data['prolif'].iloc[0]), 5)
    cell_adh = round(float(data['cell_adh'].iloc[0]), 1)
    cell_rep = data['cell_rep'].iloc[0]
    seeds = np.unique(data['seed'])

    delaunay_distance_mean = []
    delaunay_distance_std = []
    max_mot_speed = []
    r_density = []

    #### Unique simulations in the dataset
    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))

    #### Iterate over each simulation to compute statistics
    for simulation in simulations:
        df_sim = data[data['simulation'] == simulation]

        delaunay_distance = df_sim['delaunay_distance'].to_numpy()
        delaunay_distance_mean.append(np.mean(delaunay_distance))
        delaunay_distance_std.append(np.std(delaunay_distance))

        max_mot_speed.append(round(float(df_sim['max_mot_speed'].iloc[0]), 3))
        r_density.append(round(float(df_sim['ecm_density_rate'].iloc[0]), 5))

    #### Create DataFrame for heatmap
    columns = np.unique(max_mot_speed)
    index = np.flip(np.unique(r_density))
    df = pd.DataFrame(columns=columns, index=index).fillna(0.0)
    annot_df = pd.DataFrame(columns=columns, index=index).fillna('NaN')

    #### Fill DataFrame with means and standard deviations
    for mot, degr, mean, std in zip(max_mot_speed, r_density, delaunay_distance_mean, delaunay_distance_std):
        df[mot][degr] = mean
        annot = f"{mean:.2f}\n±{std:.2f}"
        annot_df[mot][degr] = annot

    annot_arr = annot_df.to_numpy()

    #### Define color palette based on ribose concentration
    if ribose == 0:
        color_rib = seaborn.color_palette('colorblind')[0]
        color_rib_dark = seaborn.color_palette('dark')[0]
    elif ribose == 50:
        color_rib = seaborn.color_palette('colorblind')[1]
        color_rib_dark = seaborn.color_palette('dark')[1]
    elif ribose == 200:
        color_rib = seaborn.color_palette('colorblind')[2]
        color_rib_dark = seaborn.color_palette('dark')[2]

    color_rib = seaborn.color_palette('colorblind')[0]
    color_rib_dark = seaborn.color_palette('dark')[0]

    #### Define colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white', color_rib, color_rib_dark])

    #### Plot heatmap
    # hmap = seaborn.heatmap(df,cmap=cmap,vmin=15, vmax=30,ax=ax)
    hmap = seaborn.heatmap(df, cmap=cmap, vmin=15, vmax=30, annot=annot_arr, annot_kws={"fontsize":10}, fmt="s")
    hmap.figure.axes[-1].yaxis.set_label_position('left')
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(label=r'Delaunay mean distance [$\mu$m]',fontsize=12)


    if contour == True:
        df = scipy.ndimage.gaussian_filter(df, 1)
        ax.contour(np.arange(.5, df.shape[1]), np.arange(.5, df.shape[0]), df, colors='white',alpha=0.5)

    #### Set axis labels and title
    plt.ylabel(r'$r_{deg,0}$ [min$^{-1}$]', color='black', fontsize=12)
    plt.yticks(color='black', fontsize=12,rotation=45,va='top')
    plt.xlabel(r'$S_0$ [$\mu$m$\cdot$min$^{-1}$]', color='black', fontsize=12)
    plt.xticks(color='black', fontsize=12)
    if title:
        # plt.title(r'$\bf{Delaunay\,mean\,distance\,at\,{%i}\,h}$' %(max(t)/60) + '\n' r'Proliferation $r_{div}$=0.00072 min$^{-1}$',color='black', fontsize=12)
        plt.title(r'Proliferation $r_{div}$=%f min$^{-1}$' %(round(prolif,3)),fontsize=12)

    plt.savefig(save_folder + f'plots/speed_vs_degr_delaunay_rib{ribose}_{simulation_name}_t{int(max(t)/60)}.png', bbox_inches="tight")
    plt.close()


def plots_speed_vs_degr_spheroid_area_growth(data, simulation_name, save_folder, title=True, contour=False):
    """
    Plot the relationship between max migration speed and degradation rate for spheroid area growth

    Parameters:
    - data: pandas DataFrame containing the simulation data.
    - simulation_name: String identifier for the simulation.
    - save_folder: Directory path to save the plot.
    - title: Boolean to decide whether to include a title on the plot
    """
    # plt.figure()
    fig, ax = plt.subplots(figsize=(7.5,5))
    seaborn.set_context("paper")
    seaborn.set_style('ticks')

    #### Collect simulation parameters from the DataFrame
    ribose = data['ribose'].iloc[0]
    t = np.unique(data['t']).astype(int)
    prolif = round(float(data['prolif'].iloc[0]), 5)
    cell_adh = round(float(data['cell_adh'].iloc[0]), 1)
    cell_rep = data['cell_rep'].iloc[0]

    spheroid_area_ratio_mean = []
    spheroid_area_ratio_std = []
    max_mot_speed = []
    r_density = []

    #### Unique simulations in the dataset
    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))

    #### Iterate over each simulation to compute statistics
    for simulation in simulations:
        df_sim = data[data['simulation'] == simulation]

        spheroid_area_init = df_sim[df_sim['t'] == min(t)]['spheroid_area'].to_numpy()
        spheroid_area_fin = df_sim[df_sim['t'] == max(t)]['spheroid_area'].to_numpy()
        spheroid_area_ratio = spheroid_area_fin / spheroid_area_init

        spheroid_area_ratio_mean.append(np.mean(spheroid_area_ratio))
        spheroid_area_ratio_std.append(np.std(spheroid_area_ratio))

        max_mot_speed.append(round(float(df_sim['max_mot_speed'].iloc[0]), 3))
        r_density.append(round(float(df_sim['ecm_density_rate'].iloc[0]), 5))

    #### Create DataFrame for heatmap
    columns = np.unique(max_mot_speed)
    index = np.flip(np.unique(r_density))
    df = pd.DataFrame(columns=columns, index=index).fillna(0.0)
    annot_df = pd.DataFrame(columns=columns, index=index).fillna('NaN')

    #### Fill DataFrame with means and standard deviations
    for mot, degr, area_mean, area_std in zip(max_mot_speed, r_density, spheroid_area_ratio_mean, spheroid_area_ratio_std):
        df[mot][degr] = area_mean
        annot = f"{area_mean:.2f}\n±{area_std:.2f}"
        annot_df[mot][degr] = annot

    annot_arr = annot_df.to_numpy()

    #### Define color palette based on ribose concentration
    if ribose == 0:
        color_rib = seaborn.color_palette('colorblind')[0]
        color_rib_dark = seaborn.color_palette('dark')[0]
    elif ribose == 50:
        color_rib = seaborn.color_palette('colorblind')[1]
        color_rib_dark = seaborn.color_palette('dark')[1]
    elif ribose == 200:
        color_rib = seaborn.color_palette('colorblind')[2]
        color_rib_dark = seaborn.color_palette('dark')[2]

    color_rib = seaborn.color_palette('colorblind')[0]
    color_rib_dark = seaborn.color_palette('dark')[0]

    #### Define colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white', color_rib, color_rib_dark])

    #### Plot heatmap
    # hmap = seaborn.heatmap(df, cmap=cmap, vmin=1, vmax=8,ax=ax)
    hmap = seaborn.heatmap(df,cmap=cmap, vmin=1, vmax=8, annot=annot_arr, annot_kws={"fontsize":10}, fmt="s")
    hmap.figure.axes[-1].yaxis.set_label_position('left')
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(label='Growth relative to t$_0$',fontsize=12)

    
    if contour == True:
        df = scipy.ndimage.gaussian_filter(df, 1)
        ax.contour(np.arange(.5, df.shape[1]), np.arange(.5, df.shape[0]), df, colors='white', alpha=0.5)

    #### Set axis labels and title
    plt.ylabel(r'$r_{deg,0}$ [min$^{-1}$]', color='black', fontsize=12)
    plt.yticks(color='black', fontsize=12,rotation=45,va='top')
    plt.xlabel(r'$S_0$ [$\mu$m$\cdot$min$^{-1}$]', color='black', fontsize=12)
    plt.xticks(color='black', fontsize=12)
    
    if title:
        # plt.title(r'$\bf{Spheroid\,growth\,relative\,to\,t_0\,at\,{%i}\,h}$' %(max(t)/60) + '\n' r'Proliferation $r_{div}$=0.00072 min$^{-1}$',color='black', fontsize=12)
        plt.title(r'Proliferation $r_{div}$=%f min$^{-1}$' %(round(prolif,3)),fontsize=12)

    plt.savefig(save_folder + f'plots/speed_vs_degr_spheroid_area_growth_rib{ribose}_{simulation_name}_t{int(max(t)/60)}.png', bbox_inches="tight")
    plt.close()


def plots_sigma_vs_delta_spheroid_area_growth(data, simulation_name, save_folder, title=True, contour=False):
    """
    Plot the relationship between sigma and delta parameters for spheroid area growth

    Parameters:
    - data: pandas DataFrame containing the simulation data.
    - simulation_name: String identifier for the simulation.
    - save_folder: Directory path to save the plot.
    - title: Boolean to decide whether to include a title on the plot
    """
    fig, ax = plt.subplots(figsize=(7.5,5))
    seaborn.set_context("paper")
    seaborn.set_style('ticks')
 
    #### Collect simulation parameters from the DataFrame
    ribose = data['ribose'].iloc[0]
    t = np.unique(data['t']).astype(int)
    prolif = round(float(data['prolif'].iloc[0]), 5)
    cell_adh = round(float(data['cell_adh'].iloc[0]), 1)
    cell_rep = data['cell_rep'].iloc[0]

    spheroid_area_ratio_mean = []
    spheroid_area_ratio_std = []
    sigma = []
    delta = []

    #### Unique simulations in the dataset
    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))

    #### Iterate over each simulation to compute statistics
    for simulation in simulations:
        df_sim = data[data['simulation'] == simulation]

        spheroid_area_init = df_sim[df_sim['t'] == min(t)]['spheroid_area'].to_numpy()
        spheroid_area_fin = df_sim[df_sim['t'] == max(t)]['spheroid_area'].to_numpy()
        spheroid_area_ratio = spheroid_area_fin / spheroid_area_init

        spheroid_area_ratio_mean.append(np.mean(spheroid_area_ratio))
        spheroid_area_ratio_std.append(np.std(spheroid_area_ratio))

        sigma.append(round(float(df_sim['sigma'].iloc[0]), 3))
        delta.append(round(float(df_sim['delta'].iloc[0]), 3))

    #### Create DataFrame for heatmap
    columns = np.unique(sigma)
    index = np.flip(np.unique(delta))
    df = pd.DataFrame(columns=columns, index=index).fillna(0.0)
    annot_df = pd.DataFrame(columns=columns, index=index).fillna('NaN')

    #### Fill DataFrame with means and standard deviations
    for a, b, area_mean, area_std in zip(sigma, delta, spheroid_area_ratio_mean, spheroid_area_ratio_std):
        df[a][b] = area_mean
        annot = f"{area_mean:.2f}\n±{area_std:.2f}"
        annot_df[a][b] = annot

    annot_arr = annot_df.to_numpy()

    #### Define color palette based on ribose concentration
    if ribose == 0:
        color_rib = seaborn.color_palette('colorblind')[0]
        color_rib_dark = seaborn.color_palette('dark')[0]
    elif ribose == 50:
        color_rib = seaborn.color_palette('colorblind')[1]
        color_rib_dark = seaborn.color_palette('dark')[1]
    elif ribose == 200:
        color_rib = seaborn.color_palette('colorblind')[2]
        color_rib_dark = seaborn.color_palette('dark')[2]

    color_rib = seaborn.color_palette('colorblind')[0]
    color_rib_dark = seaborn.color_palette('dark')[0]

    #### Define colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white', color_rib, color_rib_dark])

    #### Plot heatmap
    # hmap = seaborn.heatmap(df, cmap=cmap, vmin=1, vmax=8, ax=ax)
    hmap = seaborn.heatmap(df, cmap=cmap, vmin=1, vmax=8, annot=annot_arr, annot_kws={"fontsize":10}, fmt="s")
    hmap.figure.axes[-1].yaxis.set_label_position('left')
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=12)
    cbar.set_label(label='Growth relative to t$_0$',fontsize=12)

    if contour == True:
        df = scipy.ndimage.gaussian_filter(df, 0.7)
        ax.contour(np.arange(.5, df.shape[1]), np.arange(.5, df.shape[0]), df, colors='white')

    #### Set axis labels and title
    plt.ylabel(r"$\delta$ [mM$^{-1}$]", color='black', fontsize=12)
    plt.xlabel(r"$\sigma$ [mM$^{-1}$]", color='black', fontsize=12)
    plt.yticks(rotation=45,va='top', color='black', fontsize=12)
    plt.xticks(rotation=45,ha='right', color='black', fontsize=12)

    if title:
        # plt.title(r'$\bf{Spheroid\,growth\,relative\,to\,t_0\,at\,{%i}\,h}$' %(max(t)/60) + '\n' r'Proliferation $r_{div}$=%f min$^{-1}$' %(round(prolif,5)),fontsize=12)
        plt.title(f'Ribose {ribose} mM', fontsize=12)

    plt.savefig(save_folder + f'plots/sigma_vs_delta_spheroid_area_growth_rib{ribose}_{simulation_name}_t{int(max(t)/60)}.png', bbox_inches="tight")
    plt.close()
