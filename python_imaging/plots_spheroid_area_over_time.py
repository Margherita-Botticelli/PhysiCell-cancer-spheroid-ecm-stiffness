from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn

    
def plots_spheroid_area_over_time(data,save_folder,max_spheroid_area):
    
    #################### PLOT SPHEROID AREA OVER TIME #############################
    # print(f'Stage 0\n',flush=True)

    
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = data['t'].unique()
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    spheroid_area = data['spheroid_area'].to_numpy()

    r_anisotropy = data['anisotropy_increase_rate'].iloc[0]
    r_density = data['ecm_density_rate'].iloc[0]
    r_orientation = data['fiber_realignment_rate'].iloc[0]

    spheroid_area = np.reshape(spheroid_area, (-1,len(t))).astype(float)

    plt.figure(simulation)

    # print(f'Stage 1: data ready\n',flush=True)

    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    else: 
        color_rib = seaborn.color_palette('colorblind')[2]
    
    # for i in range(len(spheroid_area)):
    #     # Plot spheroid area over time for all random seeds
    #     plt.plot(t,spheroid_area[i],color=color_rib,alpha=0.2)

    spheroid_area_mean = np.mean(spheroid_area, axis=0)
    spheroid_area_std = np.std(spheroid_area, axis=0)

    spheroid_area_init = data[(data['t'] == 0)]['spheroid_area'].to_numpy()
    spheroid_area_fin = data[(data['t'] == 5000)]['spheroid_area'].to_numpy()

    spheroid_area_ratio = spheroid_area_fin/spheroid_area_init
    spheroid_area_ratio_mean = np.mean(spheroid_area_ratio)
    spheroid_area_ratio_std = np.std(spheroid_area_ratio)

    # print(f'Stage 2: Ready for plotting\n',flush=True)

    plt.plot(t,spheroid_area_mean,label=f'{ribose} mM',color=color_rib)
    plt.fill_between(t, spheroid_area_mean+spheroid_area_std, spheroid_area_mean-spheroid_area_std, facecolor=color_rib, alpha=0.4)
    
    # t_slope = t[20:len(t)]
    # res = stats.linregress(t_slope,spheroid_area_mean[20:len(t)])
    # y = [(i * res.slope + res.intercept) for i in t_slope]

    # plt.plot(t_slope, y, color=color_rib)

    # print(f'Stage 3: Plots ready\n',flush=True)

    ### Set axis labels
    plt.ylabel("Area [$micron^2$]",fontsize=15)
    plt.xlabel("Time [min]",fontsize=15)

    ### Set axis ticks
    plt.yticks(np.arange(0,max_spheroid_area+1000,max_spheroid_area/10),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/5),fontsize=13)

    ### Set title, overlaied plots
    # plt.suptitle(f'Spheroid area over time',fontsize=13)
    plt.title(r'$\bf{Spheroid\,area\,over\,time}$'+f'\n{prolif=}, {max_mot_speed=}, {cell_adh=}, {cell_rep=}\n{r_density=}, {r_orientation=}, {r_anisotropy=}\n{round(spheroid_area_ratio_mean, 2)}\u00B1{round(spheroid_area_ratio_std, 2)}', fontsize = 12)
    
    plt.legend(title=f'Ribose')

    # print(f'Stage 4: Plots look ready\n',flush=True)
    
    plt.savefig(save_folder + f'plots/spheroid_area_rib{ribose}_{simulation}.png', bbox_inches = "tight")

    # print(f'Stage 5: Figure spheroid_area_{simulation} saved\n',flush=True)
    # plt.close()

