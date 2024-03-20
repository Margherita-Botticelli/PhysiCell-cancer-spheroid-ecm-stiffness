from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn
import pandas as pd

def plots_adh_vs_rep(data, simulation_name, save_folder):

    #################### PLOT ADH VS REP #############################

    fig = plt.figure()

    ribose = data['ribose'].iloc[0]
    t = np.unique(data['t']).astype(float)
    prolif = data['prolif'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]

    spheroid_area_ratio_mean = []
    spheroid_area_ratio_std = []
    cell_adh = []
    cell_rep = []

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
        

        cell_adh.append(float(df_sim['cell_adh'].iloc[0]))
        cell_rep.append(float(df_sim['cell_rep'].iloc[0]))
    
    if len(np.unique(cell_adh))<len(np.unique(cell_rep)):
        size = len(np.unique(cell_adh))
    else:
        size = len(np.unique(cell_rep))

    spheroid_area_ratio = np.reshape(spheroid_area_ratio_mean, (size,-1))
    cell_rep = np.reshape(cell_rep, (size,-1))
    cell_adh = np.reshape(cell_adh, (size,-1))
    if len(np.unique(cell_adh))<len(np.unique(cell_rep)):
        ratio = [i / j for i, j in zip(cell_adh[0],cell_rep[0])]
        # ratio = [i / j for i, j in zip(cell_rep[0],cell_adh[0])]
        index = cell_adh[:,0]
    else:
        ratio = [i / j for i, j in zip(cell_adh[0],cell_rep[0])]
        index = cell_rep[:,0]
        

    #### Select color for plot depending on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    elif(ribose == 200):
        color_rib = seaborn.color_palette('colorblind')[2]
    
    
    # print(f'spheroid_area_ratio: {spheroid_area_ratio}\n', flush=True)
    
    # print(f'cell_rep: \n{cell_rep}\n', flush=True)
    # print(f'cell_adh: \n{cell_adh}\n', flush=True)
    # print(f'ratio: \n{ratio}', flush=True)

    df = pd.DataFrame(data=spheroid_area_ratio, columns=ratio, index=index)
    # df = pd.DataFrame(data=cell_rep, columns=ratio, index=cell_adh[:,0])

    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white',color_rib])

    annot_mean = ["%.2f" % number for number in spheroid_area_ratio_mean]
    annot_std = ["%.2f" % number for number in spheroid_area_ratio_std]
    annot_pm = [u'\u00B1']*len(annot_mean)
    annot_arr = np.array([a+b+c for a,b,c in zip(annot_mean,annot_pm,annot_std)])
    annot_arr = np.reshape(annot_arr, (size,-1))
    ax = seaborn.heatmap(df,cmap=cmap,vmin=1, vmax=7,annot = annot_arr, fmt="s", cbar_kws={'label': 'Growth'})

    ax.figure.axes[-1].yaxis.set_label_position('left')

    ###  Set axis labels
    if len(np.unique(cell_adh))<len(np.unique(cell_rep)):
        plt.ylabel("Adhesion strength",fontsize=15)
        plt.xlabel("Ratio repulsion/adhesion",fontsize=15)
    else:
        plt.ylabel("Repulsion strength",fontsize=15)
        plt.xlabel("Ratio adhesion/repulsion",fontsize=15)
    
    
    
    ### Set title, overlaied plots
    plt.suptitle(f'Spheroid growth relative to t$_0$',fontsize=15, y=1.0)
    plt.title(f'Prolif rate={prolif}, max mot speed={max_mot_speed}',fontsize=12)

    plt.savefig(save_folder + f'plots/adhesion_vs_repulsion_rib{ribose}_{simulation_name}.png', bbox_inches = "tight")

