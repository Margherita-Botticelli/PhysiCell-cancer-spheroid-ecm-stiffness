from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, re
import seaborn
from joblib import Parallel, delayed
from scipy import stats
from matplotlib.transforms import Affine2D
import statistics 
import pandas as pd
from skimage import io, draw, util
from tqdm import tqdm
import xml.etree.ElementTree as ET
from sklearn.cluster import DBSCAN, OPTICS


def plots_ecm_remodeling(data, simulation_name, save_folder):

    #################### PLOT ADH VS REP #############################

    fig = plt.figure()

    ribose = data['ribose'].iloc[0]
    t = data['t'].unique()
    prolif = data['prolif'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    r_anisotropy = data['anisotropy_increase_rate'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]

    spheroid_area_ratio_mean = []
    spheroid_area_ratio_std = []

    r_density = []
    r_orientation = []

    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))
    # print(f'simulations: {simulations}\n', flush=True)

    for simulation in simulations:
        spheroid_area = []

        df_sim = data[(data['simulation'] == simulation)]
        # print(f'data frame:\n {df_sim}\n', flush=True)

        # spheroid_area = 
        # print(f'spheroid_area:\n {spheroid_area}\n', flush=True)

        # spheroid_area = np.vstack(spheroid_area)
        # print(f'spheroid_area vstack:\n {spheroid_area}\n', flush=True)

        spheroid_area_init = df_sim[(df_sim['t'] == 0)]['spheroid_area'].to_numpy()
        spheroid_area_fin = df_sim[(df_sim['t'] == 5000)]['spheroid_area'].to_numpy()
        spheroid_area_ratio = spheroid_area_fin/spheroid_area_init

        spheroid_area_ratio_mean.append(np.mean(spheroid_area_ratio))
        spheroid_area_ratio_std.append(np.std(spheroid_area_ratio))
        # print(f'{spheroid_area_ratio_mean=}',flush=True)
        # print(f'{spheroid_area_ratio_std=}',flush=True)
        

        r_density.append(float(df_sim['ecm_density_rate'].iloc[0]))
        r_orientation.append(float(df_sim['fiber_realignment_rate'].iloc[0]))
    

    size = len(np.unique(r_density))

    spheroid_area_ratio = np.reshape(spheroid_area_ratio_mean, (size,-1))
    r_orientation = np.reshape(r_orientation, (size,-1))
    r_density = np.reshape(r_density, (size,-1))


    #### Select color for plot depending on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    elif(ribose == 200):
        color_rib = seaborn.color_palette('colorblind')[2]
    
    
    # print(f'spheroid_area_ratio: {spheroid_area_ratio}\n', flush=True)
    
    # print(f'r_orientation: \n{r_orientation}\n', flush=True)
    # print(f'r_density: \n{r_density}\n', flush=True)
    # print(f'ratio: \n{ratio}', flush=True)

    df = pd.DataFrame(data=spheroid_area_ratio, columns=r_density[0,:], index=r_orientation[:,0])
    # df = pd.DataFrame(data=r_orientation, columns=ratio, index=r_density[:,0])

    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white',color_rib])

    annot_mean = ["%.2f" % number for number in spheroid_area_ratio_mean]
    annot_std = ["%.2f" % number for number in spheroid_area_ratio_std]
    annot_pm = [u'\u00B1']*len(annot_mean)
    annot_arr = np.array([a+b+c for a,b,c in zip(annot_mean,annot_pm,annot_std)])
    annot_arr = np.reshape(annot_arr, (size,-1))
    ax = seaborn.heatmap(df,cmap=cmap,vmin=1, vmax=10,annot = annot_arr, fmt="s", cbar_kws={'label': 'Growth'})
    
    # ax = seaborn.heatmap(df,cmap=cmap,vmin=1, vmax=10,annot = r_orientation, cbar_kws={'label': 'Growth'})

    ax.figure.axes[-1].yaxis.set_label_position('left')

    ###  Set axis labels:
    plt.ylabel("Change in orientation rate",fontsize=12)
    plt.xlabel("Change in density rate",fontsize=12)

    
    ### Set title
    plt.title(r'$\bf{Spheroid\,growth\,relative\,to\,t_0}$'+f'\n{prolif=}, {max_mot_speed=}, {cell_adh=}, {cell_rep=}\n{r_anisotropy=}',fontsize=12)

    plt.savefig(save_folder + f'plots_ecm_remodeling_rib_{ribose}_{simulation_name}.png', bbox_inches = "tight")

