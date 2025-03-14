import sys
import numpy as np
sys.path.append(r'../python_imaging')
from pyMCDS_ECM import *

# Script REQUIRES ffmpeg to make movei!!!!!!!

######## If using on remote system, uncomment this line below to load correct matplotlib backend ################
# matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle
import imageio
import os, sys, re
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import seaborn

def print_stats(arr):
    """
    Produces relevant statistical output to screen given an array of any dimension. It flattens the in row-major style,
    the default np.flatten behavior.

    :param arr: any dimensional array, but it probably makes the most sense to be a 2-d array
    :return: Prints to termminal the array mean, quartiles, min, and max.
    """

    print("Mean: ", np.mean(arr.flatten()))
    print("Q2 quantile of arr : ", np.quantile(arr, .50))
    print("Q1 quantile of arr : ", np.quantile(arr, .25))
    print("Q3 quantile of arr : ", np.quantile(arr, .75))
    print("Min : ", arr.min())
    print("Max : ", arr.max())


def create_plot(data, snapshot, data_folder, save_name, output_plot=True, title=True):
    """
    Creates a plot as per instructions inside the function. 
    A base layer of a contour plot of either the anisotropy or the oxygen, the cells in the smulation as a scatter plot,
    and finally the ECM orientation overlaid with a quiver plot.

    Parameters
    ----------
    snapshot :
        Base name of PhysiCell output files - eg 'output00000275' --> 'output' + '%08d'
    data_folder : str
        Path to input data
    save_folder : str
        Path for image output
    output_plot : bool
        True = image file will be made. Image output is required for movie production.
    show_plot : bool
        True = plot is displayed. Expected to be false for large batches.
    Returns
    -------
    Nothing :
        Produces a png image from the input PhysiCell data.
    """
    ####################################################################################################################
    ####################################            Load data                                   ########################
    ####################################################################################################################

    #### Simulation parameters
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    r_anisotropy = data['anisotropy_increase_rate'].iloc[0]
    r_degr = data['ecm_density_rate'].iloc[0]
    # r_push = data['ecm_pushing_rate'].iloc[0]
    r_orientation = data['fiber_realignment_rate'].iloc[0]
    chemotaxis_bias = data['chemotaxis_bias'].iloc[0]
    ecm_sensitivity = data['ecm_sensitivity'].iloc[0]


    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', data_folder)

    # loads and reads ECM data
    mcds.load_ecm(snapshot + '_ECM.mat', data_folder)

    # time
    t = ( mcds.get_time() )

    # Get cell positions and attributes, microenvironment, and ECM data for plotting.

    # # Cells
    # cell_df = mcds.get_cell_df()
    # #ribose_concentration = float(cell_df['ribose_concentration'])

    #### Diffusion microenvironment
    xx, yy = mcds.get_2D_mesh()  # Mesh
    plane_oxy = mcds.get_concentrations('substrate', 0.0)  # Oxyen (used for contour plot)

    #### ECM microenvironment
    xx_ecm, yy_ecm = mcds.get_2D_ECM_mesh()  # Mesh
    dx = mcds.get_mesh_spacing()
    plane_anisotropy = mcds.get_ECM_field('anisotropy', 0.0)  # Anistropy (used for scaling and contour plot)
    # plane_anisotropy = micro # Used for contour plot
    ecm_density = mcds.get_ECM_field('density',0.0)
    ecm_density = np.reshape(ecm_density, (int(1000/dx),int(1000/dx)))
    ####################################################################################################################
    ####################################            Preprocessing                               ########################
    ####################################################################################################################

    #### Helper varialbes and functions ######

    # Number of contours (could include as a parameter)
    num_levels = 10  # 25 works well for ECM, 38 works well for oxygen

    # Make levels for contours
    levels_o2 = np.linspace(1e-14, 38, num_levels)
    # levels_ecm = np.linspace(1e-14, 1.0, num_levels)
    levels_ecm = np.linspace(0.90, 0.93, num_levels) # for the march environment - need to especially highlight small changes in anistoropy. 

    ##### Process data for plotting - weight fibers by anisotropy, mask out 0 anisotropy ECM units, get cell radii and types

    # Anisotropy strictly runs between 0 and 1. Element by element mulitplication produces weighted lengths between 0 - 1
    # for vizualization

    # scaled_ECM_x = np.multiply(mcds.data['ecm']['ECM_fields']['x_fiber_orientation'][:, :, 0], plane_anisotropy)
    # scaled_ECM_y = np.multiply(mcds.data['ecm']['ECM_fields']['y_fiber_orientation'][:, :, 0], plane_anisotropy)

    # if we want the arrows the same length instead
    ECM_x = mcds.data['ecm']['ECM_fields']['x_fiber_orientation'][:, :, 0]
    ECM_y = mcds.data['ecm']['ECM_fields']['y_fiber_orientation'][:, :, 0]

    # ECM_matrix = np.stack([ECM_x,ECM_y])
    # print(ECM_matrix.max(), ECM_matrix.min())
    # mask out zero vectors
    # mask = plane_anisotropy > 0.0001

    # # get unique cell types and radii
    # cell_df['radius'] = (cell_df['total_volume'].values * 3 / (4 * np.pi)) ** (1 / 3)

    ####################################################################################################################
    ####################################            Plotting                                    ########################
    ####################################################################################################################

    # #### Plot style
    # plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')

    seaborn.set_context("paper")
    seaborn.set_style('ticks')


    # start plot and make correct size
    # mpl.rcParams.update({'font.size': 12})
    edge = 500
    number = int(t)
    fig = plt.figure()
    ax = plt.gca()

    
    # fig, ax = plt.subplots(figsize=(12, 12))

    #### ECM density
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["white",seaborn.color_palette('colorblind')[3]])
    plt.pcolormesh(xx_ecm,yy_ecm,ecm_density[:,:],cmap=cmap,vmin=0,vmax=1)

    #### Colorbar ECM density
    cbar = plt.colorbar(shrink=0.9)
    cbar.outline.set_visible(False)
    ax.figure.axes[-1].yaxis.set_label_position('left')
    cbar.ax.tick_params(labelsize=15, colors='black') 
    cbar.set_label(label='ECM density',fontsize=15, color='black')

    # #### Oxygen
    # cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["white",seaborn.color_palette('colorblind')[4]])
    # cs = plt.contourf(xx, yy, plane_oxy, cmap=cmap, levels=levels_o2, vmin=0,vmax=38)
    # plt.colorbar(shrink=0.7,label='Chemical substrate')

    # #### Anisotropy
    # cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["white",seaborn.color_palette('colorblind')[1]])
    # cs = plt.contourf(xx_ecm, yy_ecm, plane_anisotropy, cmap="YlGnBu", levels=levels_ecm)
    # plt.pcolormesh(xx_ecm,yy_ecm,plane_anisotropy[:,:],cmap=cmap,vmin=0,vmax=1)
    # plt.colorbar(shrink=0.7,label='ECM anisotropy')

    # #### add ECM orientation vectors unscaled by anistorpy ###
    # plt.quiver(xx, yy, 20*ECM_x, 20*ECM_y,
    # pivot='middle', angles='xy', scale_units='xy', scale=1, headwidth=0,headlength=0, headaxislength=0)

    #### Add cells layer
    for j in data['ID'].tolist():
        circ = Circle((data[data['ID']==j]['position_x'].iloc[0], data[data['ID']==j]['position_y'].iloc[0]),radius=data[data['ID']==j]['radius'].iloc[0], alpha=0.7, edgecolor='black',facecolor=seaborn.color_palette('colorblind')[2])
        ax.add_artist(circ)

    #### Add scalebar 
    scalebar = AnchoredSizeBar(ax.transData,
                            100, r'100 [$\mu$m]', 'lower right', 
                            pad=0.1,
                            color='white',
                            frameon=False,
                            size_vertical=1,
                            fontproperties={'size':15})
    ax.add_artist(scalebar)


    # text_kwargs = dict(ha='center', va='top', fontsize=12, color='white')
    # plt.text(0, 450, f'Time {int(t/60)}h', **text_kwargs)

    # Labels and title (will need removed for journal - they will be added manually)
    # ax.set_xlabel(r'x [$\mu$m]',fontsize=12)
    # ax.set_ylabel(r'y [$\mu$m]',fontsize=12)
    #fig.colorbar(cs, ax=ax)
    

    # Carefully place the command to make the plot square AFTER the color bar has been added.
    ax.axis('scaled')
    # fig.tight_layout()

    plt.xticks([])#(np.arange(-edge,edge+1, step=100),fontsize=12,rotation=45)
    plt.yticks([])#(np.arange(-edge,edge+1, step=100),fontsize=12)
    plt.ylim(-edge, edge)
    plt.xlim(-edge, edge)

    if title == True:
        # plt.title(f'{prolif=}, {max_mot_speed=}\n{cell_adh=}, {cell_rep=}, {r_density=}',fontsize=10)
        plt.title(f'chemo={chemotaxis_bias}, ecm_sens={ecm_sensitivity}\nS_cm={max_mot_speed}, {r_degr=}, {r_push=}')

    #### Plot output
    if output_plot is True:
        plt.savefig(save_name,bbox_inches='tight')
    plt.close()

def create_movie(data_folder: str, save_folder: str, save_name: str):
    """
    Generates the list of files in data_folder, finds the ones with ECM data, makes plots from them, then outputs an
    ffmpeg generated movie to save_folder, naming the movie save_name.

    This function requires ffmpeg be installed at the command line.

    :param data_folder: Path to direcotry containing data
    :param save_folder: Path to save generated image and movie to
    :param save_name: Save name for movie
    :return:
    """

    # generate list of files in the directory
    files = os.listdir(data_folder)
    # files = list(filter(re.compile(r'output*ECM\.mat').search, files))

    # For all files in the directory, process only those with with 'ECM.mat' in the names. I am not sure why there is a
    # period at the beginning of the search pattern.
    for i in range(len(files)):
        if not re.search('.ECM\.mat', files[i]):
            continue

        # Sample call with meaningful variables:
        # create_plot('output00000275', data_folder='21_03_leader_follower_model_3_test/',output_plot=False, show_plot=False)
        create_plot(files[i].split('_')[0], data_folder, save_folder, output_plot=True, show_plot=False)

    # make the movie - see ffmpeg documentation for more information

    # consider not loading the unneeded data - and be sure to get rid of the unneeded fields!!!

    os.system(
        'ffmpeg -y -framerate 24 -i ' + save_folder + 'output%08d.png -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2" "' + save_name + '.mp4"')


def create_gif(data_folder: str, save_folder: str, save_name: str):
    """
    Generates the list of files in data_folder, finds the ones with ECM data, makes plots from them, then outputs 
    a gif to save_folder, naming the gif save_name.

    :param data_folder: Path to direcotry containing data
    :param save_folder: Path to save generated image and gif to
    :param save_name: Save name for gif
    :return:
    """
    images = []
    # generate list of files in the directory
    files = os.listdir(data_folder)
    # files = list(filter(re.compile(r'output*ECM\.mat').search, files))

    # For all files in the directory, process only those with with 'ECM.mat' in the names. I am not sure why there is a
    # period at the beginning of the search pattern.
    for i in range(len(files)):
        if not re.search('_ECM\.mat', files[i]):
            continue

        # Sample call with meaningful variables:
        # create_plot('output00000275', data_folder='21_03_leader_follower_model_3_test/',output_plot=False, show_plot=False)
        create_plot(files[i].split('_')[0], data_folder, save_folder, output_plot=True, show_plot=False)
        # print(save_folder + files[i].split('_')[0] + '.png')
        images.append(imageio.imread(save_folder + files[i].split('_')[0] + '.png'))
        
    # make the gif
    imageio.mimsave(save_folder + save_name + '.gif',images)


