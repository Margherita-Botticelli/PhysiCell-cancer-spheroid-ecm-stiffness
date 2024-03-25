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
from plots_adh_vs_rep import plots_adh_vs_rep
from plots_ecm_remodeling import plots_ecm_remodeling
from spheroid_area_function import spheroid_area_function
from cell_plus_environment_movie_maker import create_plot

if __name__ == '__main__':
    
    #### Figure resolution
    mpl.rcParams['figure.dpi'] = 600
    
    #### Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')

    proj = 'tests' #'ecm_density'

    #### Save folder
    save_folder = f'../results/{proj}/'  
    data_folder = f'../data/{proj}/'

    colors = seaborn.color_palette('colorblind') #+ seaborn.color_palette('dark') + seaborn.color_palette('muted') + seaborn.color_palette('bright')

    # labels = ['0mM','50mM','200mM'] #[0.0,0.0025, 0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12]
    legend_title = 'Max mot speed'# 'Adh stength'

    #### Ribose concentrations
    riboses = [0]#,50,200]

    #### Random seeds
    n_seeds = 1
    seeds = list(range(0,n_seeds))

    #### Simulations
    # table1 = [345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368]
    # table2 = [373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396]
    # table3 = [336,337,317,338,335,339,340,341,342,343,318,344,327,328,329,330,319,320,321,322,331,332,333,334]
    # table4 = list(range(401,425)) 
    # table5 = list(range(425,449)) 
    # table6 = list(range(449,473)) 

    # table7 = list(range(473,493)) 
    # table8 = list(range(493,513))
    # table9 = list(range(513,533))  
    # table10 = [513,514,515,516,517,518,519,520,573,574,575,576,521,522,523,524,577,578,579,580,581,582,583,584,585,586,587,588,525,526,527,528,589,590,591,592]

    # #### Table for varying ECM remodelling rates
    # table11 = list(range(593,618))
    # table12 = list(range(618,643))
    # table13 = list(range(643,668))
    # table14 = list(range(668,693))
    # table15 = list(range(693,718))

    # table16 = list(range(718,743))

    #### Table for ECM density only model
    table1 = list(range(0,41))

    simulations = [2]
    # simulation_name = f'{simulations[0]}_to_{simulations[-1]}'
    simulation_name = '_'.join(str(s) for s in simulations)
    
    #### Prepare parallel simulation data storage
    riboses_list = np.repeat([riboses] * len(simulations), len(seeds))
    simulations_list = np.repeat([simulations], len(seeds) * len(riboses))
    seeds_list = seeds * len(simulations) * len(riboses)
    data_folder_list = [data_folder] * len(seeds) * len(simulations) * len(riboses)

    # Parallel(n_jobs=-1)(delayed(simulation_data)(data_folder,simulation,ribose,seed) for data_folder,simulation,ribose,seed in zip(data_folder_list,simulations_list,riboses_list,seeds_list))

    simulation_data(data_folder,simulations[0],riboses[0],seeds[0])

    print('Parallel end\n', flush=True)

    #### Initiate the data frame
    df = pd.DataFrame()

    df_list = []

    #### Combine different simulations dataframes
    for sim in simulations:
        for rib in riboses:
            for seed in seeds:
                df_new = pd.read_pickle(data_folder + f'output_rib{rib}_{sim}_{seed}/dataframe_rib{rib}_{sim}_{seed}.pkl')

                # df_new = simulation_data(data_folder,sim,rib,seed) 
            
                df_list.append(df_new)

                # pd.set_option('display.max_rows', None)
                # pd.set_option('display.max_columns', None)
                # print(f'data frame new \n {df_new}', flush=True)
            

    df = pd.concat(df_list, copy=False, axis=0)
    
    ### Print dataframe by displaying max rows and columns
    # pd.set_option('display.max_rows', None)
    # pd.set_option('display.max_columns', None)

    # print('data frame \n', flush=True)
    # print(df, flush = True)

    print('Dataframe ready!\n',flush=True)


    # ########## SPHEROID AREA OVER TIME PLOT ###########
    # mpl.interactive(True)

    # max_spheroid_area = df[(df['t'] == 5000) & (df['ID'] == 0)]['spheroid_area'].max()
    # max_spheroid_area = round(max_spheroid_area,-4)+2000
    # # print(f'max_spheroid_area:{max_spheroid_area}\n',flush=True)

    # for sim in simulations:
    #     for rib in riboses:
    #         data = df[(df['simulation'] == sim) & (df['ribose'] == rib) & (df['ID'] == 0)]
    #         plots_spheroid_area_over_time(data,save_folder,max_spheroid_area)

    #     plt.close('all')
    

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
            

    # ######## ADHESION VS REPULSION HEATMAP PLOT #########
    # for rib in riboses:
    #     data = df[(df['ribose'] == rib) & (df['ID'] == 0)]
    #     plots_adh_vs_rep(data, simulation_name, save_folder)
    #     # print('Plots adh vs rep finishes\n', flush=True)
    #     plt.close('all')

    # plt.close('all')

    
    # ######## ECM REMODELING HEATMAP PLOT #########
    # for rib in riboses:
    #     data = df[(df['ribose'] == rib) & (df['ID'] == 0)]
    #     plots_ecm_remodeling(data, simulation_name, save_folder)
    #     # print('Plots adh vs rep finishes\n', flush=True)
    #     plt.close('all')

    # plt.close('all')
    

    # ######### CLUSTERING FUNCTION ##############
    # for sim in simulations:
    #     for rib in riboses:
    #         data = df[(df['simulation'] == sim) & (df['ribose'] == rib)]
    #         cluster_function(data,save_folder)
    #     # plt.close()


    ######### TIME POINT IMAGE ##############
    times = [5760]

    for sim in simulations:
        for rib in riboses:
            for t in times:
                seed = 0
                data = df[(df['simulation'] == sim) & (df['ribose'] == rib) & (df['seed'] == seed) & (df['t'] == t)]
                # spheroid_area_function(data,save_folder=save_folder,figure=True)
                # pd.set_option('display.max_columns', None)

                print('data frame \n', flush=True)
                print(data, flush = True)

                time_step = data[data['ID']==0].index.values.astype(int)[0]
                snapshot = 'output' + '{:08d}'.format(time_step)
                data_folder_sim = data_folder + f'output_rib{rib}_{sim}_{seed}/'
                save_name = save_folder + f'images/full_image_rib{rib}_{sim}_{seed}_t{int(t)}.png'

                print(f'{save_name=}\n', flush=True)
                create_plot(data, snapshot, data_folder_sim, save_name, output_plot=True, show_plot=False)


