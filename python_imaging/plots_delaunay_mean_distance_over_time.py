from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn

    
def plots_delaunay_mean_distance_over_time(data,save_folder):
    
    #################### PLOT DELAUNAY MEAN DISTANCE OVER TIME #############################
    # print(f'Stage 0\n',flush=True)

    
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = data['t'].unique() / 60
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    delaunay_distance = data['delaunay_distance'].to_numpy()

    r_anisotropy = data['anisotropy_increase_rate'].iloc[0]
    r_density = data['ecm_density_rate'].iloc[0]
    r_orientation = data['fiber_realignment_rate'].iloc[0]

    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    else: 
        color_rib = seaborn.color_palette('colorblind')[2]

    delaunay_distance = np.reshape(delaunay_distance, (-1,len(t))).astype(float)


    plt.figure(figsize=(6,3))

    # print(f'Stage 1: data ready\n',flush=True)
    
    # for i in range(len(delaunay_distance)):
    #     # Plot spheroid area over time for all random seeds
    #     plt.plot(t,delaunay_distance[i],color=color_rib,alpha=0.2)

    delaunay_distance_mean = np.mean(delaunay_distance, axis=0)
    delaunay_distance_std = np.std(delaunay_distance, axis=0)

    # print(f'Stage 2: Ready for plotting\n',flush=True)

    plt.plot(t,delaunay_distance_mean,label=f'{ribose} mM',color=color_rib)
    plt.fill_between(t, delaunay_distance_mean+delaunay_distance_std, delaunay_distance_mean-delaunay_distance_std, facecolor=color_rib, alpha=0.4)
    
    # t_slope = t[20:len(t)]
    # res = stats.linregress(t_slope,delaunay_distance_mean[20:len(t)])
    # y = [(i * res.slope + res.intercept) for i in t_slope]

    # plt.plot(t_slope, y, color=color_rib)

    # print(f'Stage 3: Plots ready\n',flush=True)

    ### Set axis labels
    plt.ylabel(f'Delaunay mean distance',fontsize=15)
    plt.xlabel("Time [h]",fontsize=15)

    ### Set axis ticks
    plt.yticks(np.arange(0,25.1,5),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/4),fontsize=13)

    ### Set title, overlaied plots
    # plt.suptitle(f'Spheroid area over time',fontsize=13)
    plt.title(r'$\bf{Delaunay\,mean\,distance}$'+f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}', fontsize = 12)
    
    plt.legend(title=f'Ribose')

    # print(f'Stage 4: Plots look ready\n',flush=True)
    
    plt.savefig(save_folder + f'plots/delaunay_mean_distance_over_time_rib{ribose}_{simulation}.png', bbox_inches = "tight")

    # print(f'Stage 5: Figure delaunay_distance_{simulation} saved\n',flush=True)
    # plt.close()

