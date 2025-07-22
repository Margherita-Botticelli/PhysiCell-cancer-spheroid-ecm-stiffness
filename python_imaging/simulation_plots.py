from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from joblib import Parallel, delayed # type: ignore
import pandas as pd
from simulation_data import *
from plots_over_time import *
from box_plots import *
from heatmaps import *
from spheroid_area_function import spheroid_area_function
from delaunay_function import *
from cell_plus_environment_movie_maker import create_plot
from skimage import io
from PIL import Image
import os
import warnings
from least_squares_function import *

#### Suppress future warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

if __name__ == '__main__':
    
    #### Set figure resolution for high-quality plots
    mpl.rcParams['figure.dpi'] = 600
    
    #### Set plotting style for consistency
    plt.style.use('seaborn-v0_8-colorblind')

    #### Specify the project to work on
    proj = 'ecm_density'  # Options: 'tests', 'ecm_fibers'
    
    #### Define folders for saving results and loading data
    save_folder = f'../results/{proj}/'
    data_folder = f'../data/{proj}/'

    #### Define simulations to run
    
    #### Non-invasive simulations
    # simulations = [240] # Set riboses = [0,50,200]
    #### Invasive simulations
    # simulations = [293] # Set riboses = [0,50,200]
    #### Invasive + MMP inhibitor simulations
    # simulations = [420] # Set riboses = [0,50,200]
    #### Heatmaps speed vs degr with r_div=0.0004, r_div=0.0008 and r_div=0.0008
    # simulations = list(range(0,240)) # Set riboses = [0]
    #### Heatmaps speed vs degr with r_div=0.00072
    simulations = list(range(240,320)) # Set riboses = [0]
    #### Heatmaps sigma vs delta with r_div=0.00072
    # simulations = list(range(320,420)) # Set riboses = [50,200]

<<<<<<< HEAD
    # simulations = list(range(0,40))   
    # simulations = list(range(0,4))
    # simulations = list(range(4,148))
    simulations = [3]

    #### Flag to determine if existing data should be replaced
    replace = True
    
    #### List of ribose concentrations to test
    # riboses = [0] #[50,200] # [50,200] # Options: 0, 50, 200

    #### List of fibre orientations to test
    orientations = ['random','radial', 'tangential'] # ['tangential'] # Options:'random', 'radial', 'tangential'
=======
    # simulations = list(range(0,3))
    # simulations = [240,293]

    #### Flag to determine if existing data should be replaced
    replacing = False
    
    #### List of ribose concentrations to test
    riboses = [0] #[50,200] # [50,200] # Options: 0, 50, 200
>>>>>>> parent of 7c561b3 (Update simulation_plots.py)

    #### Number of random seeds for simulations
    n_seeds = 1
    seeds = list(range(0, n_seeds))

    #### Flags for different types of plots
<<<<<<< HEAD
    title = True # Title on plots
    least_squares = False
    plots_over_time = True 
=======
    title = False # Title on plots
    least_squares = True
    plots_over_time = False 
>>>>>>> parent of 7c561b3 (Update simulation_plots.py)
    box_plots = False
    contour = False # Contour lines on heatmaps
    heatmaps_speed_vs_degr = False
    heatmaps_sigma_vs_delta = False
    heatmap_time_points = [72*60]  # Example: [48*60], [96*60]
    time_point_images = False # 'all' # True # Options: 'all' if you want all time points
    times = [72*60] #[96*60] # Time points to consider for images, Options: [24*60, 48*60, 72*60, 96*60]
    video = False
    statistics = False

    #### Define a name for the simulation set based on its indices
    simulation_name = '_'.join(str(s) for s in simulations)
    
    ########## PROCESSING SIMULATION DATA ###########
    #### Prepare lists for parallel simulation data
    riboses_list = np.repeat([riboses] * len(simulations), len(seeds))
    simulations_list = np.repeat([simulations], len(seeds) * len(riboses))
    seeds_list = seeds * len(simulations) * len(riboses)
    data_folder_list = [data_folder] * len(seeds) * len(simulations) * len(riboses)
    replace_list = [replacing] * len(seeds) * len(simulations) * len(riboses)

    # simulations_list = np.repeat([simulations], len(seeds))
    # seeds_list = seeds * len(simulations) 
    # data_folder_list = [data_folder] * len(seeds) * len(simulations)
    # replace_list = [replacing] * len(seeds) * len(simulations) 

    ### Run simulations in parallel
    # Parallel(n_jobs=10)(delayed(simulation_data)(data_folder, simulation, ribose, seed, replace)
    #                         for data_folder, simulation, ribose, seed, replace in zip(data_folder_list, simulations_list, riboses_list, seeds_list, replace_list))

    # Parallel(n_jobs=10)(delayed(simulation_data)(data_folder, simulation, seed, replace) 
    #                         for data_folder, simulation, seed, replace in zip(data_folder_list, simulations_list, seeds_list, replace_list))

    #### Notify when parallel processing is complete
    # print('Parallel end\n', flush=True)

    #### Initialize an empty DataFrame and list to store simulation data
    df = pd.DataFrame()
    df_list = []

    #### Combine data from different simulations into one DataFrame
    for sim in simulations:
        for rib in riboses:
            for seed in seeds:
                df_new = pd.read_pickle(data_folder + f'output_rib{rib}_{sim}_{seed}/dataframe_rib{rib}_{sim}_{seed}.pkl')
                # df_new = pd.read_pickle(data_folder + f'output_{sim}_{seed}/dataframe_{sim}_{seed}.pkl')
                df_list.append(df_new)

    #### Concatenate all dataframes into a single dataframe
    df = pd.concat(df_list, copy=False, axis=0)

    print('Dataframe ready!\n', flush=True)
    
    if least_squares:
        data = df[df['ID'] == 0]
        # least_squares_function_approximated(data, save_folder)
        least_squares_function_cell_count(df, save_folder)

    ########## PLOTS OVER TIME ###########
    if plots_over_time:
        for sim in simulations:
            for rib in riboses:
                data = df[(df['simulation'] == sim) & (df['ribose'] == rib)]
                print(f'plots over time {sim=}', flush=True)

<<<<<<< HEAD
                    data = df[(df['simulation'] == sim) & (df['orientation'] == orientation)]

                    plots_spheroid_area_growth_over_time(data, save_folder, title=title)
                    plots_delaunay_mean_distance_over_time(data, save_folder, title=title)
                    # plots_cell_count_over_time(data, save_folder, title=title)
                    plots_invasion_over_time(data, save_folder, n, title=title)
                plt.close('all')
=======
                # data = df[(df['simulation'] == sim)]
                # plots_spheroid_area_growth_over_time(data, save_folder, title=title)
                # plots_delaunay_mean_distance_over_time(data, save_folder, title=title)
                plots_cell_count_over_time(data, save_folder, title=title)
                # plots_invasion_over_time(data, save_folder, title=title)
                # plt.close('all')
            plt.close('all')
>>>>>>> parent of 7c561b3 (Update simulation_plots.py)

        print('plots over time done!', flush=True)


    ########## BOX PLOTS ###########
    if box_plots:
        if sim == 420:
            box_plot_spheroid_area_growth(df, save_folder, title=title)
            plt.close('all')
        else:
            for sim in simulations:
                data = df[(df['simulation'] == sim)]
                # box_plot_spheroid_area_growth(data, save_folder, title=title)
                box_plot_delaunay_mean_distance(data, save_folder, title=title)
                box_plot_cell_count(data, save_folder, title=title)
            plt.close('all')

        print('Box plots done!', flush=True)


    ######## HEATMAP PLOTS #########

    #### Initialize empty DataFrames for specific plots
    data_spheroid_area_growth = pd.DataFrame()
    data_cluster = pd.DataFrame()
    data_delaunay = pd.DataFrame()
    data_cell_count = pd.DataFrame()
    
    for rib in riboses:
        for time_point in heatmap_time_points:

            #### Check if heatmaps for speed vs degradation should be generated
            if heatmaps_speed_vs_degr:

                #### Filter data to loop through unique proliferation rates
                data = df[(df['ribose'] == rib) & (df['ID'] == 0) & ((df['t'] == 0) | (df['t'] == time_point))]
                for prolif in np.unique(data['prolif']).astype(float):

                    #### Further filter data by the current proliferation rate
                    data_spheroid_area_growth = df[(df['ribose'] == rib) & (df['ID'] == 0) & ((df['t'] == 0) | (df['t'] == time_point)) & (df['prolif'] == prolif)]
                    
                    #### Get first simulation ID to name the plot
                    simulations_heatmap = np.unique(data_spheroid_area_growth['simulation']).astype(int)
                    simulation_name_plot = simulations_heatmap[0]
                    
                    #### Generate the heatmap for spheroid area growth
                    plots_speed_vs_degr_spheroid_area_growth(data_spheroid_area_growth, simulation_name_plot, save_folder, title=title, contour=contour)
                    print('plots_speed_vs_degr_spheroid_area_growth done!', flush=True)

                #### Filter data to loop through unique proliferation rates
                data = df[(df['ribose'] == rib) & (df['ID'] == 0) & (df['t'] == time_point)]
                for prolif in np.unique(data['prolif']).astype(float):

                    #### Further filter data by the current proliferation rate
                    data_delaunay = data[(data['prolif'] == prolif)]
                    
                    #### Get first simulation ID to name the plot
                    simulations_heatmap = np.unique(data_delaunay['simulation']).astype(int)
                    simulation_name_plot = simulations_heatmap[0]
                    
                    #### Generate the heatmap for Delaunay mean distance
                    plots_speed_vs_degr_delaunay(data_delaunay, simulation_name_plot, save_folder, title=title, contour=contour)
                    print('plots_speed_vs_degr_delaunay done!', flush=True)

            #### Check if heatmaps for sigma vs delta should be generated
            if heatmaps_sigma_vs_delta:
                #### Filter data
                data_spheroid_area_growth = df[(df['ribose'] == rib) & (df['ID'] == 0) & ((df['t'] == 0) | (df['t'] == time_point))]
                
                #### Get first simulation ID to name the plot
                simulations_heatmap = np.unique(data_spheroid_area_growth['simulation']).astype(int)
                simulation_name_plot = simulations_heatmap[0]
                
                #### Generate the heatmap for spheroid area growth
                plots_sigma_vs_delta_spheroid_area_growth(data_spheroid_area_growth, simulation_name_plot, save_folder, title=title, contour=contour)
                print('plots_sigma_vs_delta_spheroid_area_growth done!', flush=True)

    #### Close all plot figures to free up memory
    plt.close('all')


    ######### TIME POINT IMAGE ##############
    #### Check if time point images should be generated
    if(time_point_images or video):
        for sim in simulations:
            for seed in seeds:
                for rib in riboses:
                    if time_point_images:
                        if time_point_images == 'all':
                            times = np.unique(df[(df['seed'] == seed) & (df['ID'] == 0)]['t']).astype(int)
                            # print('dataframe', flush=True)
                        for t in times:

                            #### Filter data
                            seed = 0
                            data = df[(df['simulation'] == sim) & (df['ribose'] == rib) & (df['seed'] == seed) & (df['t'] == t)]
                            # data = df[(df['simulation'] == sim) & (df['seed'] == seed) & (df['t'] == t)]

                            # print(data, flush=True)
                            #### Get time point to find snapshot
                            time_step = data[data['ID'] == 0].index.values.astype(int)[0]
                            snapshot = 'output' + '{:08d}'.format(int(time_step))
                            data_folder_sim = data_folder + f'output_rib{rib}_{sim}_{seed}/'
                            save_name = save_folder + f'images/full_image_rib{rib}_{sim}_{seed}_t{int(t):04}.png'
                            # data_folder_sim = data_folder + f'output_{sim}_{seed}/'
                            # save_name = save_folder + f'images/full_image_{sim}_{seed}_t{int(t):04}.png'

                            print(f'{sim=}, {t=}', flush=True)
                            #### Generate images
                            create_plot(data, snapshot, data_folder_sim, save_name, output_plot=True, title=title) 
                            plt.close('all')

                
                #### Check if video should be generated
                if video:
                    # video_name = save_folder + f'animations/video_rib{rib}_{sim}_{seed}.mp4'
                    video_name = save_folder + f'animations/video_{sim}_{seed}.mp4'

                    #### Find generated images
                    # images = save_folder + f'images/full_image_rib{rib}_{sim}_{seed}_t*.png'
                    images = save_folder + f'images/full_image_{sim}_{seed}_t*.png'

                    #### Generate video
                    os.system(f'ffmpeg -y -framerate 10 -pattern_type glob -i \'{images}\' -c:v libx264 -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" {video_name}')

                    print('Video ready!', flush=True)

    ######### STATISTICS IMAGES ##########
    if statistics:
        for sim in simulations:
            data = df[(df['simulation'] == sim) & (df['seed'] == 0) & (df['t'] == 96*60)]
            delaunay_distance_function(data, save_folder=save_folder, figure=True)
            spheroid_area_function(data, save_folder=save_folder, figure=True)