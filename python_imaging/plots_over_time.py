from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn

    
def plots_spheroid_growth_over_time(data,save_folder,title=True):
    
    #################### PLOT SPHEROID AREA OVER TIME #############################
    # print(f'Stage 0\n',flush=True)

    data = data[(data['ID'] == 0)]

    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = data['t'].unique() / 60
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    spheroid_area = data['spheroid_area'].to_numpy()

    r_anisotropy = data['anisotropy_increase_rate'].iloc[0]
    r_density = data['ecm_density_rate'].iloc[0]
    r_orientation = data['fiber_realignment_rate'].iloc[0]

    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    else: 
        color_rib = seaborn.color_palette('colorblind')[2]

    spheroid_area = np.reshape(spheroid_area, (-1,len(t))).astype(float)

    spheroid_area_ratio = spheroid_area[:,0]
    spheroid_area_ratio = np.reshape(spheroid_area_ratio,(-1,1))

    spheroid_area_ratio = spheroid_area / spheroid_area_ratio

    plt.figure(figsize=(5,4),num=simulation)

    # print(f'Stage 1: data ready\n',flush=True)
    
    # for i in range(len(spheroid_area)):
    #     # Plot spheroid area over time for all random seeds
    #     plt.plot(t,spheroid_area[i],color=color_rib,alpha=0.2)

    spheroid_area_mean = np.mean(spheroid_area_ratio, axis=0)
    spheroid_area_std = np.std(spheroid_area_ratio, axis=0)

    # print(f'Stage 2: Ready for plotting\n',flush=True)

    plt.plot(t,spheroid_area_mean,label=f'{ribose} mM',color=color_rib)
    plt.fill_between(t, spheroid_area_mean+spheroid_area_std, spheroid_area_mean-spheroid_area_std, facecolor=color_rib, alpha=0.4)
    
    # t_slope = t[20:len(t)]
    # res = stats.linregress(t_slope,spheroid_area_mean[20:len(t)])
    # y = [(i * res.slope + res.intercept) for i in t_slope]

    # plt.plot(t_slope, y, color=color_rib)

    # print(f'Stage 3: Plots ready\n',flush=True)

    ### Set axis labels
    plt.ylabel(f'Growth relative to t$_0$',fontsize=18)
    plt.xlabel("Time [h]",fontsize=18)

    ### Set axis ticks
    plt.yticks(np.arange(0,8.1,1),fontsize=15)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/4),fontsize=15)

    ### Set title, overlaied plots
    if title == True:
        # plt.suptitle(f'Spheroid area over time',fontsize=15)
        plt.title(r'$\bf{Spheroid\,growth\,relative\,to\,t_0}$'+f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}', fontsize = 12)
    
    plt.legend(title=f'Ribose',fontsize=15,title_fontsize=15)

    # print(f'Stage 4: Plots look ready\n',flush=True)
    
    plt.savefig(save_folder + f'plots/spheroid_growth_over_time_{simulation}.png', bbox_inches = "tight")

    # print(f'Stage 5: Figure spheroid_area_{simulation} saved\n',flush=True)
    # plt.close()

   
def plots_delaunay_mean_distance_over_time(data,save_folder, title=True):
    
    #################### PLOT DELAUNAY MEAN DISTANCE OVER TIME #############################
    # print(f'Stage 0\n',flush=True)
    data = data[(data['ID'] == 0)]
    
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


    plt.figure(figsize=(5,4),num=simulation+1)

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
    plt.ylabel(f'Delaunay mean distance',fontsize=18)
    plt.xlabel("Time [h]",fontsize=18)

    ### Set axis ticks
    plt.yticks(np.arange(0,25.1,5),fontsize=15)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/4),fontsize=15)

    ### Set title, overlaied plots
    if title == True:
        # plt.suptitle(f'Spheroid area over time',fontsize=15)
        plt.title(r'$\bf{Delaunay\,mean\,distance}$'+f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}', fontsize = 12)
        
    plt.legend(title=f'Ribose',fontsize=15,title_fontsize=15)
    

    # print(f'Stage 4: Plots look ready\n',flush=True)
    
    plt.savefig(save_folder + f'plots/delaunay_mean_distance_over_time_{simulation}.png', bbox_inches = "tight")

    # print(f'Stage 5: Figure delaunay_distance_{simulation} saved\n',flush=True)
    # plt.close()


   
def plots_cell_number_over_time(data,save_folder,title=True):
    
    #################### PLOT CELL NUMBER OVER TIME #############################
    # print(f'Stage 0\n',flush=True)

    
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = data['t'].unique() / 60
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    seeds = data['seed'].unique()

    r_anisotropy = data['anisotropy_increase_rate'].iloc[0]
    r_density = data['ecm_density_rate'].iloc[0]
    r_orientation = data['fiber_realignment_rate'].iloc[0]

    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    else: 
        color_rib = seaborn.color_palette('colorblind')[2]


    cell_number = []
    for seed in seeds:
        for timepoint in data['t'].unique(): 
        # print(f'{timepoint=}', flush=True)
            # print(f'{seed=}', flush=True)
            df_seed = data[(data['t'] == timepoint) & (data['seed'] == seed)]
            cell_number_t = df_seed['ID'].iloc[-1]
            cell_number.append(cell_number_t+1)

    cell_number = np.reshape(cell_number, (-1,len(t))).astype(int)
    # print(f'{cell_number=}', flush=True)



    cell_number_mean = np.mean(cell_number,axis=0)
    cell_number_std = np.std(cell_number,axis=0)

    plt.figure(figsize=(5,4),num=simulation+2)

    plt.plot(t,cell_number_mean,label=f'{ribose} mM',color=color_rib)
    plt.fill_between(t, cell_number_mean+cell_number_std, cell_number_mean-cell_number_std, facecolor=color_rib, alpha=0.4)
    
    ### Set axis labels
    plt.ylabel(f'Cell number',fontsize=18)
    plt.xlabel("Time [h]",fontsize=18)

    ### Set axis ticks
    plt.yticks(np.arange(0,1200,200),fontsize=15)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/4),fontsize=15)

    ### Set title, overlaied plots
    if title == True:
        # plt.suptitle(f'Spheroid area over time',fontsize=15)
        plt.title(r'$\bf{Cell\,number}$'+f'\n{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}', fontsize = 12)
        
    plt.legend(title=f'Ribose',fontsize=15,title_fontsize=15)
    


    # print(f'Stage 4: Plots look ready\n',flush=True)
    
    plt.savefig(save_folder + f'plots/cell_number_over_time_{simulation}.png', bbox_inches = "tight")

    # print(f'Stage 5: Figure delaunay_distance_{simulation} saved\n',flush=True)
    # plt.close()

