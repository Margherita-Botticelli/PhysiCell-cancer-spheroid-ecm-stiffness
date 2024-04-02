from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn
from joblib import Parallel, delayed
import pandas as pd
from simulation_data import *
from plots_spheroid_area_over_time import plots_spheroid_area_over_time
from cluster_function import cluster_function
from plots_adh_vs_rep import *
from plots_ecm_remodeling import plots_ecm_remodeling
from spheroid_area_function import spheroid_area_function
from cell_plus_environment_movie_maker import create_plot
from cell_velocities import cell_velocities
from plots_spheroid_growth_over_time import plots_spheroid_growth_over_time
from skimage import io
from PIL import Image


if __name__ == '__main__':
    
    #### Figure resolution
    mpl.rcParams['figure.dpi'] = 600
    
    #### Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')

    #### Table for ECM density only model
    table1 = list(range(0,25))
    table2 = list(range(25,50))
    table3 = list(range(50,75))

    table4 = list(range(75,100))
    table5 = list(range(100,125))
    table6 = list(range(125,150))

    table7 = list(range(150,175))
    table8 = list(range(175,200))
    table9 = list(range(200,225))

    #### Choose simulations
    simulations_multi = [table1,table2,table3,table4,table5,table6,table7,table8,table9]
    
    #### Select project
    proj = 'ecm_density' # 'tests' # 

    #### Save folder and data folder
    save_folder = f'../results/{proj}/'  
    data_folder = f'../data/{proj}/'

    #### Ribose concentrations
    riboses = [0]#,50,200]

    #### Random seeds
    n_seeds = 5
    seeds = list(range(0,n_seeds))

    #### Replacing dataframes True or False
    replacing = False

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
        # print(df, flush = True)

        print('Dataframe ready!\n',flush=True)


        # ########## SPHEROID GROWTH OVER TIME PLOT ###########

        # for sim in simulations:
        #     for rib in riboses:
        #         data = df[(df['simulation'] == sim) & (df['ribose'] == rib) & (df['ID'] == 0)]
        #         plots_spheroid_growth_over_time(data,save_folder)

        #     plt.close('all')
        

        ######## ADHESION VS REPULSION HEATMAP PLOT #########
        for rib in riboses:
            data_spheroid_growth = df[(df['ribose'] == rib) & (df['ID'] == 0)]
            plots_adh_vs_rep_spheroid_growth(data_spheroid_growth, simulation_name, save_folder)

            data_cluster = df[(df['ribose'] == rib) & (df['t'] == 5760)]
            plots_adh_vs_rep_clusters(data_cluster, simulation_name, save_folder)
            # print('Plots adh vs rep finishes\n', flush=True)
            plt.close('all')

        plt.close('all')

        
        # ######## ECM REMODELING HEATMAP PLOT #########
        # for rib in riboses:
        #     data = df[(df['ribose'] == rib) & (df['ID'] == 0)]
        #     plots_ecm_remodeling(data, simulation_name, save_folder)
        #     # print('Plots adh vs rep finishes\n', flush=True)
        #     plt.close('all')

        # plt.close('all')
        

        # ######### TIME POINT IMAGE ##############
        # times = np.unique(df[(df['seed'] == 0) & (df['ID'] == 0)]['t']).astype(int) # [1440,2880,5760]
        
        # for sim in simulations:
        #     for rib in riboses:
        #         images = []
        #         for t in times:
        #             seed = 0
        #             data = df[(df['simulation'] == sim) & (df['ribose'] == rib) & (df['seed'] == seed) & (df['t'] == t)]
        #             # spheroid_area_function(data,save_folder=save_folder,figure=True)
        #             # pd.set_option('display.max_columns', None)

        #             # print('data frame \n', flush=True)
        #             # print(data, flush = True)

        #             time_step = data[data['ID']==0].index.values.astype(int)[0]
        #             snapshot = 'output' + '{:08d}'.format(time_step)
        #             data_folder_sim = data_folder + f'output_rib{rib}_{sim}_{seed}/'
        #             save_name = save_folder + f'images/full_image_rib{rib}_{sim}_{seed}_t{int(t)}.png'

        #             print(f'{save_name=}\n', flush=True)
        #             create_plot(data, snapshot, data_folder_sim, save_name, output_plot=True, show_plot=False)

        #             images.append(save_name)
                
        #         imgs = []
        #         for i in images:
        #             print(i)
        #             imgs.append(io.imread(i))
        #         animation = [Image.fromarray(i) for i in imgs]
        #         animation[0].save(save_folder + f'animations/video_rib{rib}_{sim}_0.gif', save_all=True, append_images=animation[1:], duration=50, loop=0)
        #         del animation


        # ######### CELL VELOCITIES ##############
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

        
        # ########## SPHEROID AREA GROWTH FOR DIFFERENT PARAMETER SETS ############
        # mpl.interactive(True)
        # for sim in simulations:
        #     for rib in riboses:
        #         data = df[(df['simulation'] == sim) & (df['ribose'] == rib) & (df['ID'] == 0)]

        #         ##### plot for different adh/rep ratios
        #         cell_rep = data['cell_rep'].values.tolist()
        #         cell_adh = data['cell_adh'].values.tolist()
        #         adh_rep_ratio = float(cell_adh[0])/float(cell_rep[0])
        #         plots_spheroid_area_growth(data,int(adh_rep_ratio*100),save_folder,simulation_name)

        #         ##### plot for different max motility speeds
        #         # plots_spheroid_area_growth(data,sum(simulations),save_folder,simulation_name)

        # plt.close('all')


        # ########## SPHEROID AREA OVER TIME PLOT ###########
        # mpl.interactive(True)

        # max_spheroid_area = df[(df['t'] == 5760) & (df['ID'] == 0)]['spheroid_area'].max()
        # max_spheroid_area = round(max_spheroid_area,-4)+2000
        # # print(f'max_spheroid_area:{max_spheroid_area}\n',flush=True)

        # for sim in simulations:
        #     for rib in riboses:
        #         data = df[(df['simulation'] == sim) & (df['ribose'] == rib) & (df['ID'] == 0)]
        #         plots_spheroid_area_over_time(data,save_folder,max_spheroid_area)

        #     plt.close('all')
        