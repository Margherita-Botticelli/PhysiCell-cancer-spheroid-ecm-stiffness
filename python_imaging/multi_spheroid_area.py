from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys, os, re, math
import imageio
from matplotlib.patches import Circle
from matplotlib import cm
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

def spheroid_area_data(snapshot, data_folder, save_folder):
    print(flush=True)
    #################################  Load data  ############################
    
    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', data_folder)

    # loads and reads ECM data
    mcds.load_ecm(snapshot + '_ECM.mat', data_folder)

    # time
    t = ( mcds.get_time() )

    # Cells
    cell_df = mcds.get_cell_df()

    for j in cell_df.index:
        voxel_centre_j = voxel_centre(cell_df.loc[j, 'position_x'], cell_df.loc[j, 'position_y'])
        voxel_centres.append(voxel_centre_j)

    # print(voxel_centres,flush=True)

    # Remove doubles
    voxel_centres = list(set(voxel_centres))

    # Compute area
    voxel_side = 20
    voxel_number = len(list)
    spheroid_area = voxel_number * voxel_side**2

    number_cells = len(cell_df.index)
    # print(voxel_centres,flush=True)

    # x_val = [x[0] for x in voxel_centres]
    # y_val = [x[1] for x in voxel_centres]

    # plt.plot(x_val,y_val,'ro')

    return t, spheroid_area, number_cells


def spheroid_area_plot(ribose,simulation):

    # Save folder
    save_folder = '../results/'  #results' + simulation + '/'     

    
    # Initiate figure, single plot
    #plt.figure()
    spheroid_area = []
    number_cells = []
    results = []

    if ribose == 0:
        color = seaborn.color_palette('colorblind')[0]
    elif ribose == 50:
        color = seaborn.color_palette('colorblind')[1]
    elif ribose == 200:
        color = seaborn.color_palette('colorblind')[3]

    # Number of simulations with different random seed
    seeds = list(range(6))

    simulations = [simulation] * len(seeds)
    # print(f"simulations {simulations}")
    riboses = [ribose] * len(seeds)
    # print(f"riboses {riboses}")

    results = list(zip(*Parallel(n_jobs=6)(delayed(spheroid_area_data)(r,s,seed) for r,s,seed in zip(riboses,simulations,seeds))))

    # print(results)
    t = results[0][0]
    spheroid_area = results[1]
    number_cells = results[2]
    # print(t)
    # print(spheroid_area)
    
    plt.figure()
   
    for i in range(len(spheroid_area)):
        # Plot average migration speeds over time
        plt.plot(t,spheroid_area[i],color=color,alpha=0.2)

    spheroid_area = np.mean(spheroid_area, axis=0)
    zorder = 100 + ribose
    plt.plot(t,spheroid_area,label=f'{ribose}mM',color=color,zorder=zorder)

    #  Set axis labels
    plt.ylabel("Area $[micron^2]$",fontsize=15)
    plt.xlabel("Time [min]",fontsize=15)
    
    print({t[-1]})
    # Set axis ticks
    plt.yticks(np.arange(0,151,20),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/5),fontsize=13)
    #plt.xticks(np.arange(0, 2000+1,200),fontsize=13)
    
    # Set title, single plot
    #plt.title('Average cells total speeds, ribose {ribose}mM'.format(ribose=ribose))
    
    # Set title, overlaied plots
    plt.title('Area convered by spheroid')
    
    plt.legend()
    
    plt.savefig(save_folder + 'spheroid_area_' + str(simulation) + '.png', bbox_inches = "tight")

    plt.close()

    plt.figure()
   
    for i in range(len(number_cells)):
        # Plot average migration speeds over time
        plt.plot(t,number_cells[i],color=color,alpha=0.2)

    number_cells = np.mean(number_cells, axis=0)
    zorder = 100 + ribose
    plt.plot(t,number_cells,label=f'{ribose}mM',color=color,zorder=zorder)

    #  Set axis labels
    plt.ylabel("Number of cells",fontsize=15)
    plt.xlabel("Time [min]",fontsize=15)
    
    print({t[-1]})
    # Set axis ticks
    plt.yticks(np.arange(0,151,20),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/5),fontsize=13)
    #plt.xticks(np.arange(0, 2000+1,200),fontsize=13)
    
    # Set title, single plot
    #plt.title('Average cells total speeds, ribose {ribose}mM'.format(ribose=ribose))
    
    # Set title, overlaied plots
    plt.title('Total number of cells')
    
    plt.legend()
    
    plt.savefig(save_folder + 'number_cells_' + str(simulation) + '.png', bbox_inches = "tight")

    plt.close()

    #plt.show()


if __name__ == '__main__':
    
    # Figure resolution
    mpl.rcParams['figure.dpi'] = 300
    
    # Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')
    
    # Initiate figure
    plt.figure()
    
    # Ribose concentrations
    ribose = [0,50,200]
    # ribose = ['50']
    
    simulation = 10
    # simulations = [simulation] * len(ribose)

    for ribose in ribose:
        # print(list(zip(ribose, simulations)),flush=True)
        spheroid_area_plot(ribose,simulation)

    # Parallel(n_jobs=3)(delayed(max_dist_plot)(r,s) for r,s in zip(ribose,simulations))

    # for r in ribose:
    #     # Max distance reached w.r.t time
    #     max_dist_plot(r,simulation)
        

    plt.close()
    # fig = plt.gcf()
    # fig.patch.set_alpha(0.0)