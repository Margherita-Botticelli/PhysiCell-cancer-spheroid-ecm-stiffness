from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, re
import seaborn
from joblib import Parallel, delayed

def voxel_centre(x,y):
    # x += 10
    x /= 20
    # sign_x = np.sign(x)
    x = np.floor(x)  * 20 + 10
    # j = int(x)
    # y = -y
    # y += 10
    y /= 20
    # sign_y = np.sign(y)
    y = np.floor(y) * 20 + 10
    # i = int(y)

    return x,y

def snapshot_data(snapshot,data_folder):
    
    ### Function for average cells speeds over time 
    
    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', data_folder)

    # time
    t = ( mcds.get_time() )

    # Cells
    cell_df = mcds.get_cell_df()
    
    # Get each cell's speed
    dist = []
    voxel_centres = []

    for j in cell_df.index:
        position_x = cell_df.loc[j, 'position_x']
        position_y = cell_df.loc[j, 'position_y']
        
        dist.append(np.sqrt(position_x**2 + position_y**2))

        voxel_centre_j = voxel_centre(cell_df.loc[j, 'position_x'], cell_df.loc[j, 'position_y'])
        voxel_centres.append(voxel_centre_j)

        speed_single = cell_df.loc[j, 'total_speed']

    # print(voxel_centres,flush=True)

    # Remove doubles
    voxel_centres = list(set(voxel_centres))

    # Compute area
    voxel_side = 20
    voxels_number = len(voxel_centres)
    spheroid_area = voxels_number * voxel_side**2

    number_cells = len(cell_df.index)
    
    max_dist = max(dist)

    return t, max_dist, spheroid_area, number_cells, speed_single

def simulation_data(ribose,simulation,seed):
 
    print(f"Simulation {simulation}_{seed} with ribose {ribose}\n", flush=True)

    # Data folder
    data_folder = '../ribose_' + str(ribose) + '/output/output_' + str(simulation)  + '_' + str(seed) + '/'

    t = []
    max_dist = []
    spheroid_area = []
    number_cells = []
    speed_single = []
    

    files = os.listdir(data_folder)
    files.sort()
    # print(files)
    
    for i in range(len(files)):
        if not re.search('_ECM\.mat', files[i]):
            continue
        
        # Call function to get time and average total speeds values
        t_i, max_dist_i, spheroid_area_i, number_cells_i, speed_single_i = snapshot_data(files[i].split('_')[0], data_folder)
        
        t.append(t_i)
        
        max_dist.append(max_dist_i)
        spheroid_area.append(spheroid_area_i)
        number_cells.append(number_cells_i)
        speed_single.append(speed_single_i)

    # max_dist_mean.append(max_dist)
     
    return t, max_dist, spheroid_area, number_cells, speed_single

def plots(ribose,simulation,n_seeds):

    # Save folder
    save_folder = '../results/'  #results' + simulation + '/'     

    
    # Initiate figure, single plot
    #plt.figure()

    results = []
    max_dist_mean = []
    spheroid_area = []
    number_cells = []
    speed_single = []


    if ribose == 0:
        color = seaborn.color_palette('colorblind')[0]
    elif ribose == 50:
        color = seaborn.color_palette('colorblind')[1]
    elif ribose == 200:
        color = seaborn.color_palette('colorblind')[3]

    
    # Number of simulations with different random seed
    seeds = list(range(n_seeds))
    seeds.remove(2)

    simulations = [simulation] * len(seeds)
    # print(f"simulations {simulations}")
    riboses = [ribose] * len(seeds)
    # print(f"riboses {riboses}")

    results = list(zip(*Parallel(n_jobs=len(seeds))(delayed(simulation_data)(r,s,seed) for r,s,seed in zip(riboses,simulations,seeds))))

    print(f'', flush=True)

    # print(results)
    t = results[0][0]
    max_dist_mean = results[1]
    spheroid_area = results[2]
    number_cells = results[3]
    speed_single = results[4]

    # print(t)
    # print(max_dist_mean)
    


    #################### PLOT MAX DISTANCE #############################
    plt.figure(1)
    
    for i in range(len(max_dist_mean)):
        # Plot average migration speeds over time
        plt.plot(t,max_dist_mean[i],color=color,alpha=0.2)

    max_dist_mean = np.mean(max_dist_mean, axis=0)
    zorder = 100 + ribose
    plt.plot(t,max_dist_mean,label=f'{ribose}mM',color=color,zorder=zorder)

    #  Set axis labels
    plt.ylabel("Distance [micron]",fontsize=15)
    plt.xlabel("Time [min]",fontsize=15)
    
    # print({t[-1]})
    # Set axis ticks
    plt.yticks(np.arange(0,501,20),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/5),fontsize=13)

    # Set title, single plot
    #plt.title('Average cells total speeds, ribose {ribose}mM'.format(ribose=ribose))
    
    # Set title, overlaied plots
    plt.title('Max distance from origin')
    
    plt.legend()
    
    plt.savefig(save_folder + 'multi_max_dist_' + str(simulation) + '.png', bbox_inches = "tight")

    # plt.close()

    #################### PLOT SPHEROID AREA #############################
    plt.figure(2)
    
    for i in range(len(spheroid_area)):
        # Plot average migration speeds over time
        plt.plot(t,spheroid_area[i],color=color,alpha=0.2)

    spheroid_area = np.mean(spheroid_area, axis=0)
    zorder = 100 + ribose
    plt.plot(t,spheroid_area,label=f'{ribose}mM',color=color,zorder=zorder)

    #  Set axis labels
    plt.ylabel("Area [$micron^2$]",fontsize=15)
    plt.xlabel("Time [min]",fontsize=15)
    
    # print({t[-1]})
    # Set axis ticks
    # plt.yticks(np.arange(0,100,20),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/5),fontsize=13)

    # Set title, single plot
    #plt.title('Average cells total speeds, ribose {ribose}mM'.format(ribose=ribose))
    
    # Set title, overlaied plots
    plt.title('Area covered by the spheroid')
    
    plt.legend()
    
    plt.savefig(save_folder + 'spheroid_area_' + str(simulation) + '.png', bbox_inches = "tight")

    
    # #################### PLOT SPEED SINGLE CELL #############################

    # plt.figure(3)
    
    # for i in range(len(speed_single)):
    #     # Plot average migration speeds over time
    #     plt.plot(t,speed_single[i],color=color,alpha=0.2)

    # speed_single = np.mean(speed_single, axis=0)
    # zorder = 100 + ribose
    # plt.plot(t,speed_single,label=f'{ribose}mM',color=color,zorder=zorder)

    # #  Set axis labels
    # plt.ylabel("Speed",fontsize=15)
    # plt.xlabel("Time [min]",fontsize=15)
    
    # # print({t[-1]})
    # # Set axis ticks
    # # plt.yticks(np.arange(0,300,10),fontsize=13)
    # plt.xticks(np.arange(0, t[-1]+1,t[-1]/5),fontsize=13)

    # # Set title, single plot
    # #plt.title('Average cells total speeds, ribose {ribose}mM'.format(ribose=ribose))
    
    # # Set title, overlaied plots
    # plt.title('Speed of single cell')
    
    # plt.legend()
    
    # plt.savefig(save_folder + 'speed_single_' + str(simulation) + '.png', bbox_inches = "tight")
    
    #################### PLOT NUMBER OF CELLS #############################

    plt.figure(4)
    
    for i in range(len(number_cells)):
        # Plot average migration speeds over time
        plt.plot(t,number_cells[i],color=color,alpha=0.2)

    number_cells = np.mean(number_cells, axis=0)
    zorder = 100 + ribose
    plt.plot(t,number_cells,label=f'{ribose}mM',color=color,zorder=zorder)

    #  Set axis labels
    plt.ylabel("Number of cells",fontsize=15)
    plt.xlabel("Time [min]",fontsize=15)
    
    # print({t[-1]})
    # Set axis ticks
    # plt.yticks(np.arange(0,300,10),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/5),fontsize=13)

    # Set title, single plot
    #plt.title('Average cells total speeds, ribose {ribose}mM'.format(ribose=ribose))
    
    # Set title, overlaied plots
    plt.title('Total number of cells')
    
    plt.legend()

    plt.savefig(save_folder + 'number_cells_' + str(simulation) + '.png', bbox_inches = "tight")

    # plt.close()
    mpl.interactive(True)

if __name__ == '__main__':
    
    # Figure resolution
    mpl.rcParams['figure.dpi'] = 300
    
    # Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')

    
    # Ribose concentrations
    riboses = [0,50,200]
    # ribose = ['50']
    
    simulations = [87] #range(81,85) 

    n_seeds = 6

    for simulation in simulations:
        print(f'==== SIMULATION {simulation} STARTED ====\n', flush=True)
                
        # Initiate figure
        plt.figure()

        for ribose in riboses:

            print(f'==== RIBOSE {ribose} STARTED ====\n', flush=True)

            plots(ribose,simulation,n_seeds)

            print(f'==== RIBOSE {ribose} DONE ====\n', flush=True)

        plt.close('all')

        print(f'==== SIMULATION {simulation} ENDED ====\n', flush=True)

        # fig = plt.gcf()
        # fig.patch.set_alpha(0.0)