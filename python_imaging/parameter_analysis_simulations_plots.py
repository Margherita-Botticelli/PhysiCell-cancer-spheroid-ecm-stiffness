from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn
from joblib import Parallel, delayed
import pandas as pd
from simulation_data import *
from cluster_function import cluster_function
from heatmaps import *
from plots_ecm_remodeling import plots_ecm_remodeling
from spheroid_area_function import spheroid_area_function
from delaunay_function import *
from cell_plus_environment_movie_maker import create_plot
from cell_velocities import cell_velocities
from plots_over_time import *
from skimage import io
from PIL import Image
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


if __name__ == '__main__':
    
    #### Figure resolution
    mpl.rcParams['figure.dpi'] = 600
    
    #### Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')

    #### Select project
    proj = 'ecm_density' # 'tests' # 'ecm_fibres' # 

    #### Save folder and data folder
    save_folder = f'../results/{proj}/'  
    data_folder = f'../data/{proj}/'

    #### Choose simulations
    # simulations_multi = table2
    simulations_multi = [[240]] # invasive
    # simulations_multi = [[291]] # non-invasive
    # simulations_multi = [list(range(0,240))]
    # simulations_multi = [list(range(242,245))]
    # simulations_multi = [list(range(241,290))]
  
    #### Replacing dataframes True or False
    replacing = False
    
    #### Ribose concentrations
    riboses = [0]#,50,200] #[0]#,50,200]

    #### Random seeds, indicate number of seeds
    n_seeds = 1
    seeds = list(range(0,n_seeds))

    #### Plots
    plots_over_time = False
    title = False
    heatmaps = False
    time_point_images = False
    times = 'all' # [96*60] # [24*60,96*60] # [24*60] # [24*60,48*60,72*60,96*60] # [24*60,96*60] # 
    video = False
    statistics = False


    for simulations in simulations_multi:
        #### Define simulation name for figures
        simulation_name = '_'.join(str(s) for s in simulations)
        
        #### Prepare parallel simulations
        riboses_list = np.repeat([riboses] * len(simulations), len(seeds))
        simulations_list = np.repeat([simulations], len(seeds) * len(riboses))
        seeds_list = seeds * len(simulations) * len(riboses)
        data_folder_list = [data_folder] * len(seeds) * len(simulations) * len(riboses)
        replace_list = [replacing] * len(seeds) * len(simulations) * len(riboses)

        #### Run parallel simulations
        Parallel(n_jobs=10)(delayed(simulation_data)(data_folder,simulation,ribose,seed,replace) for data_folder,simulation,ribose,seed,replace in zip(data_folder_list,simulations_list,riboses_list,seeds_list,replace_list))

        # #### Run single simulation
        # simulation_data(data_folder,simulations[0],riboses[0],seeds[0])

        print('Parallel end\n', flush=True)

        #### Initiate the data frame
        df = pd.DataFrame()
        df_list = []

        #### Combine different simulations dataframes
        for sim in simulations:
            for rib in riboses:
                for seed in seeds:
                    df_new = pd.read_pickle(data_folder + f'output_rib{rib}_{sim}_{seed}/dataframe_rib{rib}_{sim}_{seed}.pkl')
                    
                    df_list.append(df_new)

                    # pd.set_option('display.max_rows', None)
                    # pd.set_option('display.max_columns', None)
                    # print(f'data frame new \n {df_new}', flush=True)

        df = pd.concat(df_list, copy=False, axis=0)
        
        # #### Print dataframe by displaying max rows and columns
        # pd.set_option('display.max_rows', None)
        # pd.set_option('display.max_columns', None)

        # print('data frame \n', flush=True)
        # data = df[(df['t']==0) & (df['ID']==0)]
        # print(data[['simulation','prolif', 'max_mot_speed', 'cell_adh', 'cell_rep', 'ecm_density_rate']], flush = True)
        # # print(df.columns, flush = True)
        # print(df, flush = True)

        print('Dataframe ready!\n',flush=True)

        ########## PLOTS OVER TIME ###########
        if plots_over_time==True:
            for sim in simulations:
                for rib in riboses:
                    data = df[(df['simulation'] == sim) & (df['ribose'] == rib)]
                    plots_spheroid_growth_over_time(data,save_folder,title)
                    plots_delaunay_mean_distance_over_time(data,save_folder,title)
                    plots_cell_number_over_time(data,save_folder,title)

                plt.close('all')

        data_spheroid_growth = pd.DataFrame()
        data_cluster = pd.DataFrame()
        data_delaunay = pd.DataFrame()
        data_cell_number = pd.DataFrame()

        # heatmap_num = 0
        
        ######## HEATMAP PLOTS #########
        if heatmaps==True:
            timepoints = [96*60]
            for rib in riboses:
                for timepoint in timepoints:
                        # data = df[(df['ribose'] == rib) & (df['ID'] == 0) & ((df['t'] == 0) | (df['t'] == timepoint))]
                        # for prolif in np.unique(data['prolif']).astype(float):
                        #     data_spheroid_growth =df[(df['ribose'] == rib) & (df['ID'] == 0) & ((df['t'] == 0) | (df['t'] == timepoint)) & (df['prolif'] == prolif)]
                        #     simulations_heatmap = np.unique(data_spheroid_growth['simulation']).astype(int)
                        #     simulation_name_plot = simulations_heatmap[0]
                        #     plots_mot_vs_degr_spheroid_growth(data_spheroid_growth, simulation_name_plot, save_folder, title)
                        #     print('plots_mot_vs_degr_spheroid_growth done!',flush=True)
                            

                        # data = df[(df['ribose'] == rib) & (df['ID'] == 0) & (df['t'] == timepoint)]
                        # for prolif in np.unique(data['prolif']).astype(float):
                        #     data_delaunay = data[(data['prolif'] == prolif)]
                        #     simulations_heatmap = np.unique(data_delaunay['simulation']).astype(int)
                        #     # simulation_name_plot = '_'.join(str(s) for s in simulations)
                        #     simulation_name_plot = simulations_heatmap[0]
                            
                        #     plots_mot_vs_degr_delaunay(data_delaunay, simulation_name_plot, save_folder, title)
                        #     print('plots_mot_vs_degr_delaunay done!',flush=True)                   

                        
                        data_spheroid_growth = df[(df['ribose'] == rib) & (df['ID'] == 0) & ((df['t'] == 0) | (df['t'] == timepoint))]
                        simulations_heatmap = np.unique(data_spheroid_growth['simulation']).astype(int)
                        # simulation_name_plot = '_'.join(str(s) for s in simulations)
                        simulation_name_plot = simulations_heatmap[0]
                        plots_alpha_vs_beta_spheroid_growth(data_spheroid_growth, simulation_name_plot, save_folder, title)
                        print('plots_alpha_vs_beta_spheroid_growth done!',flush=True)
                            
                        

                    # data_spheroid_growth = df[(df['ribose'] == rib) & (df['ID'] == 0) & ((df['t'] == 0) | (df['t'] == timepoint))]
                    # plots_adh_vs_rep_spheroid_growth(data_spheroid_growth, simulation_name, save_folder)
                    # print('plots_adh_vs_rep_spheroid_growth done!',flush=True)

                    # # data_cluster = df[(df['ribose'] == rib) & (df['t'] == 5760)]
                    # # plots_adh_vs_rep_clusters(data_cluster, simulation_name, save_folder)
                    # # print('plots_adh_vs_rep_clusters done!',flush=True)

                    # data_delaunay = df[(df['ribose'] == rib) & (df['t'] == timepoint)]
                    # plots_adh_vs_rep_delaunay(data_delaunay, simulation_name, save_folder)
                    # print('plots_adh_vs_rep_delaunay done!',flush=True)

                    # data_cell_number = df[(df['ribose'] == rib) & (df['t'] == timepoint)]
                    # plots_adh_vs_rep_cell_number(data_cell_number, simulation_name, save_folder)
                    # print('plots_adh_vs_rep_cell_number done!',flush=True)

            plt.close('all')

        
        # ######## ECM REMODELING HEATMAP PLOT #########
        # for rib in riboses:
        #     data = df[(df['ribose'] == rib) & (df['ID'] == 0)]
        #     plots_ecm_remodeling(data, simulation_name, save_folder)
        #     # print('Plots adh vs rep finishes\n', flush=True)
        #     plt.close('all')

        # plt.close('all')

        ######### TIME POINT IMAGE ##############
        if times == 'all':
            times = np.unique(df[(df['seed'] == 0) & (df['ID'] == 0)]['t']).astype(int)
        
        for sim in simulations:
            print(f'{sim=}',flush=True)
            for rib in riboses:
                if time_point_images==True:
                    for t in times:

                        # print(sim,rib,t, flush = True)
                        seed = 0
                        data = (df[(df['simulation'] == sim) & (df['ribose'] == rib) & (df['seed'] == seed) & (df['t'] == t)])

                        # spheroid_area_function(data,save_folder=save_folder,figure=True)
                        # pd.set_option('display.max_columns', None)

                        # print('data frame \n', flush=True)
                        # print(data, flush = True)

                        time_step = data[data['ID']==0].index.values.astype(int)[0]
                        snapshot = 'output' + '{:08d}'.format(int(time_step))
                        data_folder_sim = data_folder + f'output_rib{rib}_{sim}_{seed}/'
                        save_name = save_folder + f'images/full_image_rib{rib}_{sim}_{seed}_t{int(t):04}.png'

                        create_plot(data, snapshot, data_folder_sim, save_name, output_plot=True, title=title) 
                    
                if video==True:
                    seed = 0
                    video_name = save_folder + f'animations/video_rib{rib}_{sim}_{seed}.mp4'
                    images = save_folder + f'images/full_image_rib{rib}_{sim}_{seed}_t*.png'

                    os.system(f'ffmpeg -y -framerate 10 -pattern_type glob -i \'{images}\' -c:v libx264 -pix_fmt yuv420p -vf pad=\"width=ceil(iw/2)*2:height=ceil(ih/2)*2\" {video_name}')



        # ########### CELL VELOCITIES ##############
        # for sim in simulations:
        #     for rib in riboses:
        #         data = df[(df['simulation'] == sim) & (df['ribose'] == rib)]
        #         cell_velocities(data,save_folder)
        #     # plt.close()

        
        # ######### CLUSTERING FUNCTION ##############
        # for sim in simulations:
        #     for rib in riboses:
        #         for seed in seeds:
        #             times = np.unique(df[(df['seed'] == 0) & (df['ID'] == 0)]['t']).astype(int) # [1440,2880,5760]
        #             for t in times:
        #                 data = df[(df['simulation'] == sim) & (df['ribose'] == rib) & (df['seed'] == seed) & (df['t'] == t)]
        #                 cluster_function(data,save_folder,figure=False)
        #     # plt.close()

        if statistics==True:
            data = df[df['t'] == 96*60]
            print(data,flush=True)
            delaunay_distance_function(data,save_folder=save_folder,figure=True)
            spheroid_area_function(data,save_folder=save_folder,figure=True)