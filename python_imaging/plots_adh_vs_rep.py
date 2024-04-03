from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn
import pandas as pd
from cluster_function import cluster_function
from delaunay_function import delaunay_distance_function

def plots_adh_vs_rep_spheroid_growth(data, simulation_name, save_folder):

    #################### PLOT ADH VS REP FOR SPHEROID GROWTH #############################

    plt.figure()

    ribose = data['ribose'].iloc[0]
    t = np.unique(data['t']).astype(int)
    print(f'{t=}',flush=True)
    prolif = round(float(data['prolif'].iloc[0]),5)
    max_mot_speed = data['max_mot_speed'].iloc[0]

    spheroid_area_ratio_mean = []
    spheroid_area_ratio_std = []
    cell_adh = []
    cell_rep = []

    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))
    # print(f'simulations: {simulations}\n', flush=True)

    for simulation in simulations:
        df_sim = data[(data['simulation'] == simulation)]
        # print(f'data frame:\n {df_sim}\n', flush=True)

        spheroid_area_init = df_sim[(df_sim['t'] == min(t))]['spheroid_area'].to_numpy()
        spheroid_area_fin = df_sim[(df_sim['t'] == max(t))]['spheroid_area'].to_numpy()
        spheroid_area_ratio = spheroid_area_fin/spheroid_area_init

        spheroid_area_ratio_mean.append(np.mean(spheroid_area_ratio))
        spheroid_area_ratio_std.append(np.std(spheroid_area_ratio))
        # print(f'{spheroid_area_ratio_mean=}',flush=True)
        # print(f'{spheroid_area_ratio_std=}',flush=True)

        cell_adh.append(float(df_sim['cell_adh'].iloc[0]))
        cell_rep.append(float(df_sim['cell_rep'].iloc[0]))
    
    columns = np.unique(cell_adh) 
    index = np.flip(np.unique(cell_rep))

    df = pd.DataFrame(columns=columns, index=index)
    df = df.fillna(0.0)

    annot_df = pd.DataFrame(columns=columns, index=index)
    annot_df = annot_df.fillna('NaN')
    # print(annot_df,flush=True)
    for adh,rep,area_mean,area_std in zip(cell_adh,cell_rep,spheroid_area_ratio_mean,spheroid_area_ratio_std):
        df[adh][rep] = area_mean

        annot_mean = "%.2f" % area_mean
        annot_std = "%.2f" % area_std
        annot = annot_mean + '\n' + u'\u00B1' + annot_std
        annot_df[adh][rep] = annot
    
    # print(annot_df,flush=True)

    annot_arr = annot_df.to_numpy()
    # print(annot_arr,flush=True)

    #### Select color for plot depending on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
        color_rib_dark = seaborn.color_palette('dark')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
        color_rib_dark = seaborn.color_palette('dark')[1]
    elif(ribose == 200):
        color_rib = seaborn.color_palette('colorblind')[2]
        color_rib_dark = seaborn.color_palette('dark')[2]

    # df = pd.DataFrame(data=spheroid_area_ratio, columns=cell_adh[0,:], index=cell_rep[:,0])
    # df = pd.DataFrame(data=cell_rep, columns=ratio, index=cell_adh[:,0])

    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white',color_rib,color_rib_dark])

    ax = seaborn.heatmap(df,cmap=cmap,vmin=0, vmax=11,annot = annot_arr, fmt="s", cbar_kws={'label': 'Growth'})
    # ax = seaborn.heatmap(df,cmap=cmap,vmin=0, vmax=7,annot=True,fmt='.2f', cbar_kws={'label': 'Growth'})

    ax.figure.axes[-1].yaxis.set_label_position('left')

    ###  Set axis labels
    plt.ylabel("Repulsion strength",fontsize=15)
    plt.xlabel("Adhesion strength",fontsize=15)

    ### Set title, overlaied plots
    plt.title(r'$\bf{Spheroid\,growth\,relative\,to\,t_0}$'+
              f'\nProlif rate={prolif}, max mot speed={max_mot_speed}' + 
              f'\nt={max(t)/60}h', fontsize = 12)
    
    plt.savefig(save_folder + f'plots/adh_vs_rep_spheroid_growth_rib{ribose}_{simulation_name}_t{int(max(t)/60)}.png', bbox_inches = "tight")
    plt.close()


def plots_adh_vs_rep_clusters(data, simulation_name, save_folder):

    #################### PLOT ADH VS REP FOR MAIN CLUSTER PERCENTAGE #############################
    
    plt.figure()

    #### Collect the data from the dataframe
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = np.unique(data['t']).astype(int)
    # seed = data['seed'].iloc[0]
    prolif = round(float(data['prolif'].iloc[0]),5)
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    
    seeds = np.unique(data['seed'])
    # print(f'{seeds=}',flush=True)

    main_cluster_cell_percentages_mean = []
    main_cluster_cell_percentages_std = []
    cell_adh = []
    cell_rep = []

    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))
    # print(f'simulations: {simulations}\n', flush=True)

    for simulation in simulations:
        df_sim = data[(data['simulation'] == simulation)]
        main_cluster_cell_percentages = []
        for seed in seeds:
            df_seed = df_sim[(df_sim['seed'] == seed)]
            if seed == 0:
                figure = True
            else: figure = False
            main_cluster_cell_percentages.append(cluster_function(df_seed,save_folder,figure=figure))
 
        main_cluster_cell_percentages_mean.append(np.mean(main_cluster_cell_percentages))
        # print(f'{main_cluster_cell_percentages_mean=}',flush=True)
        main_cluster_cell_percentages_std.append(np.std(main_cluster_cell_percentages))

        cell_adh.append(float(df_sim['cell_adh'].iloc[0]))
        cell_rep.append(float(df_sim['cell_rep'].iloc[0]))

    columns = np.unique(cell_adh) 
    index = np.flip(np.unique(cell_rep))

    df = pd.DataFrame(columns=columns, index=index)
    df = df.fillna(0.0)

    annot_df = pd.DataFrame(columns=columns, index=index)
    annot_df = annot_df.fillna('NaN')
    # print(annot_df,flush=True)
    for adh,rep,cluster_mean,cluster_std in zip(cell_adh,cell_rep,main_cluster_cell_percentages_mean,main_cluster_cell_percentages_std):
        df[adh][rep] = cluster_mean

        annot_mean = "%.2f" % cluster_mean
        annot_std = "%.2f" % cluster_std
        annot = annot_mean + '\n' + u'\u00B1' + annot_std
        annot_df[adh][rep] = annot
    
    annot_arr = annot_df.to_numpy()

    #### Select color for plot depending on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
        color_rib_dark = seaborn.color_palette('dark')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
        color_rib_dark = seaborn.color_palette('dark')[1]
    elif(ribose == 200):
        color_rib = seaborn.color_palette('colorblind')[2]
        color_rib_dark = seaborn.color_palette('dark')[2]

    # df = pd.DataFrame(data=spheroid_area_ratio, columns=cell_adh[0,:], index=cell_rep[:,0])
    # df = pd.DataFrame(data=cell_rep, columns=ratio, index=cell_adh[:,0])

    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white',color_rib,color_rib_dark])

    # print(df,flush=True)


    ax = seaborn.heatmap(df,cmap=cmap,vmin=0, vmax=100,annot = annot_arr, fmt="s", cbar_kws={'label': 'Growth'})
    # ax = seaborn.heatmap(df,cmap=cmap,vmin=0, vmax=7,annot=True,fmt='.2f', cbar_kws={'label': 'Growth'})

    ax.figure.axes[-1].yaxis.set_label_position('left')

    ###  Set axis labels
    plt.ylabel("Repulsion strength",fontsize=15)
    plt.xlabel("Adhesion strength",fontsize=15)
    # if len(np.unique(cell_adh))<len(np.unique(cell_rep)):
    #     plt.ylabel("Adhesion strength",fontsize=15)
    #     plt.xlabel("Ratio repulsion/adhesion",fontsize=15)
    # else:
    #     plt.ylabel("Repulsion strength",fontsize=15)
    #     plt.xlabel("Ratio adhesion/repulsion",fontsize=15)
    
    
    ### Set title, overlaied plots
    plt.title(r'$\bf{Percentage\,of\,cells\,in\,main\,cluster}$'+
              f'\nProlif rate={prolif}, max mot speed={max_mot_speed}' + 
              f'\nt={max(t)/60}h', fontsize = 12)
    

    plt.savefig(save_folder + f'plots/adh_vs_rep_cluster_rib{ribose}_{simulation_name}_t{int(max(t)/60)}.png', bbox_inches = "tight")

    plt.close()


def plots_adh_vs_rep_delaunay(data, simulation_name, save_folder):

    #################### PLOT ADH VS REP FOR MAIN DELAUNAY #############################
    
    plt.figure()

    #### Collect the data from the dataframe
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = np.unique(data['t']).astype(int)
    # seed = data['seed'].iloc[0]
    prolif = round(float(data['prolif'].iloc[0]),5)
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    
    seeds = np.unique(data['seed'])
    # print(f'{seeds=}',flush=True)

    delaunay_distance_mean = []
    delaunay_distance_std = []
    cell_adh = []
    cell_rep = []

    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))
    # print(f'simulations: {simulations}\n', flush=True)

    for simulation in simulations:
        df_sim = data[(data['simulation'] == simulation)]
        delaunay_distance = []
        for seed in seeds:
            df_seed = df_sim[(df_sim['seed'] == seed)]
            if seed == 0:
                figure = True
            else: figure = False
            delaunay_distance.append(delaunay_distance_function(df_seed,save_folder,figure=figure))
 
        delaunay_distance_mean.append(np.mean(delaunay_distance))
        # print(f'{main_cluster_cell_percentages_mean=}',flush=True)
        delaunay_distance_std.append(np.std(delaunay_distance))

        cell_adh.append(float(df_sim['cell_adh'].iloc[0]))
        cell_rep.append(float(df_sim['cell_rep'].iloc[0]))

    columns = np.unique(cell_adh) 
    index = np.flip(np.unique(cell_rep))

    df = pd.DataFrame(columns=columns, index=index)
    df = df.fillna(0.0)

    annot_df = pd.DataFrame(columns=columns, index=index)
    annot_df = annot_df.fillna('NaN')
    # print(annot_df,flush=True)
    for adh,rep,cluster_mean,cluster_std in zip(cell_adh,cell_rep,delaunay_distance_mean,delaunay_distance_std):
        df[adh][rep] = cluster_mean

        annot_mean = "%.2f" % cluster_mean
        annot_std = "%.2f" % cluster_std
        annot = annot_mean + '\n' + u'\u00B1' + annot_std
        annot_df[adh][rep] = annot
    
    annot_arr = annot_df.to_numpy()

    #### Select color for plot depending on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
        color_rib_dark = seaborn.color_palette('dark')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
        color_rib_dark = seaborn.color_palette('dark')[1]
    elif(ribose == 200):
        color_rib = seaborn.color_palette('colorblind')[2]
        color_rib_dark = seaborn.color_palette('dark')[2]

    # df = pd.DataFrame(data=spheroid_area_ratio, columns=cell_adh[0,:], index=cell_rep[:,0])
    # df = pd.DataFrame(data=cell_rep, columns=ratio, index=cell_adh[:,0])

    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white',color_rib,color_rib_dark])

    # print(df,flush=True)


    ax = seaborn.heatmap(df,cmap=cmap,vmin=0, vmax=35,annot = annot_arr, fmt="s", cbar_kws={'label': 'Growth'})
    # ax = seaborn.heatmap(df,cmap=cmap,vmin=0, vmax=7,annot=True,fmt='.2f', cbar_kws={'label': 'Growth'})

    ax.figure.axes[-1].yaxis.set_label_position('left')

    ###  Set axis labels
    plt.ylabel("Repulsion strength",fontsize=15)
    plt.xlabel("Adhesion strength",fontsize=15)
    # if len(np.unique(cell_adh))<len(np.unique(cell_rep)):
    #     plt.ylabel("Adhesion strength",fontsize=15)
    #     plt.xlabel("Ratio repulsion/adhesion",fontsize=15)
    # else:
    #     plt.ylabel("Repulsion strength",fontsize=15)
    #     plt.xlabel("Ratio adhesion/repulsion",fontsize=15)
    
    
    ### Set title, overlaied plots
    plt.title(r'$\bf{Mean\,Delaunay\,distance}$' +
              f'\nProlif rate={prolif}, max mot speed={max_mot_speed}' + 
              f'\nt={max(t)/60}h', fontsize = 12)


    plt.savefig(save_folder + f'plots/adh_vs_rep_delaunay_rib{ribose}_{simulation_name}_t{int(max(t)/60)}.png', bbox_inches = "tight")
    plt.close()




def plots_adh_vs_rep_cell_number(data, simulation_name, save_folder):

    #################### PLOT ADH VS REP FOR NUMBER F CELLS #############################
    
    plt.figure()

    #### Collect the data from the dataframe
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    t = np.unique(data['t']).astype(int)
    # seed = data['seed'].iloc[0]
    prolif = round(float(data['prolif'].iloc[0]),5)
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    
    seeds = np.unique(data['seed'])
    # print(f'{seeds=}',flush=True)

    cell_number_mean = []
    cell_number_std = []
    cell_adh = []
    cell_rep = []

    simulations = list(dict.fromkeys(data['simulation'].values.tolist()))
    # print(f'simulations: {simulations}\n', flush=True)

    for simulation in simulations:
        df_sim = data[(data['simulation'] == simulation)]
        cell_number = []
        for seed in seeds:
            df_seed = df_sim[(df_sim['seed'] == seed)]

            cell_number.append(df_seed['ID'].iloc[-1]+1)

        cell_number_mean.append(np.mean(cell_number))
        # print(f'{main_cluster_cell_percentages_mean=}',flush=True)
        cell_number_std.append(np.std(cell_number))

        cell_adh.append(float(df_sim['cell_adh'].iloc[0]))
        cell_rep.append(float(df_sim['cell_rep'].iloc[0]))

    columns = np.unique(cell_adh) 
    index = np.flip(np.unique(cell_rep))

    df = pd.DataFrame(columns=columns, index=index)
    df = df.fillna(0.0)

    annot_df = pd.DataFrame(columns=columns, index=index)
    annot_df = annot_df.fillna('NaN')
    # print(annot_df,flush=True)
    for adh,rep,cluster_mean,cluster_std in zip(cell_adh,cell_rep,cell_number_mean,cell_number_std):
        df[adh][rep] = cluster_mean

        annot_mean = "%.2f" % cluster_mean
        annot_std = "%.2f" % cluster_std
        annot = annot_mean + '\n' + u'\u00B1' + annot_std
        annot_df[adh][rep] = annot
    
    annot_arr = annot_df.to_numpy()

    #### Select color for plot depending on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
        color_rib_dark = seaborn.color_palette('dark')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
        color_rib_dark = seaborn.color_palette('dark')[1]
    elif(ribose == 200):
        color_rib = seaborn.color_palette('colorblind')[2]
        color_rib_dark = seaborn.color_palette('dark')[2]

    # df = pd.DataFrame(data=spheroid_area_ratio, columns=cell_adh[0,:], index=cell_rep[:,0])
    # df = pd.DataFrame(data=cell_rep, columns=ratio, index=cell_adh[:,0])

    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white',color_rib,color_rib_dark])

    # print(df,flush=True)


    ax = seaborn.heatmap(df,cmap=cmap,vmin=0, vmax=1500,annot = annot_arr, fmt="s", cbar_kws={'label': 'Number of cells'})
    # ax = seaborn.heatmap(df,cmap=cmap,vmin=0, vmax=7,annot=True,fmt='.2f', cbar_kws={'label': 'Growth'})

    ax.figure.axes[-1].yaxis.set_label_position('left')

    ###  Set axis labels
    plt.ylabel("Repulsion strength",fontsize=15)
    plt.xlabel("Adhesion strength",fontsize=15)
    
    ### Set title, overlaied plots
    plt.title(r'$\bf{Number\,of\,cells}$'+
              f'\nProlif rate={prolif}, max mot speed={max_mot_speed}' + 
              f'\nt={max(t)/60}h', fontsize = 12)
    
    plt.savefig(save_folder + f'plots/adh_vs_rep_cell_number_rib{ribose}_{simulation_name}_t{int(max(t)/60)}.png', bbox_inches = "tight")
    plt.close()


