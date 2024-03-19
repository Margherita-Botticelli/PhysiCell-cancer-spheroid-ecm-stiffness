from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys, os, re, math
import imageio
from matplotlib.patches import Circle
from matplotlib.patches import Rectangle
from matplotlib import cm
import seaborn
from joblib import Parallel, delayed

def voxel_centre(x,y,size):
    # x += 10
    x /= size
    # sign_x = np.sign(x)
    x = np.floor(x)  * size # + size/2
    # j = int(x)
    # y = -y
    # y += 10
    y /= size
    # sign_y = np.sign(y)
    y = np.floor(y) * size #+ size/2
    # i = int(y)

    return x,y

def cell_spheroid_data(snapshot, data_folder, save_folder):
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

    # load cell and microenvironment data
    mcds_init = pyMCDS('initial.xml', data_folder)
    cell_df_init = mcds_init.get_cell_df()

    # Cells
    max_mot_speed = cell_df_init.loc[0,'migration_speed']
    prolif = cell_df_init.loc[0,'current_cycle_phase_exit_rate']
    # rib = cell_df.loc[0,'ribose_concentration']
    cell_adh = cell_df_init.loc[0,'cell_cell_adhesion_strength']
    cell_rep = cell_df_init.loc[0,'cell_cell_repulsion_strength']


    # # Get ribose concetration
    # ribose_concentration = float(cell_df.loc[cell_df.index[0],'ribose_concentration'])

    #### Diffusion microenvironment
    xx, yy = mcds.get_2D_mesh()  # Mesh

    #### ECM microenvironment
    xx_ecm, yy_ecm = mcds.get_2D_ECM_mesh()  # Mesh
    # plane_anisotropy = mcds.get_ECM_field('anisotropy', 0.0)  # Anistropy (used for scaling and contour plot)
    ecm_density = mcds.get_ECM_field('density',0.0)
    ecm_density = np.reshape(ecm_density, (50,50))

    #############################  Preprocessing  ########################

    #### Helper varialbes and functions ######

    # Number of contours (could include as a parameter)
    num_levels = 25  # 25 works well for ECM, 38 works well for oxygen

    # Make levels for contours
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

    # mask out zero vectors
    # mask = plane_anisotropy > 0.0001

    # get unique cell types and radii
    cell_df['radius'] = (cell_df['total_volume'].values * 3 / (4 * np.pi)) ** (1 / 3)

    ####################################  Plotting    ########################

    ##### start plot and make correct size
    mpl.rcParams.update({'font.size': 22})
    edge = 200
    fig, ax = plt.subplots(figsize=(12, 12))
    plt.ylim(-edge, edge)
    plt.xlim(-edge, edge)
    
    #### ECM density colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["#fff5f0",'#B86B56'])
    # cmap = 'copper_r'

    #### Plot ECM density mech
    #plt.contourf(X,Y,ecm_density[:,:],levels=levels,cmap=cmap)
    plt.pcolormesh(xx_ecm,yy_ecm,ecm_density[:,:],cmap=cmap,vmin=0,vmax=1)
    
    #### ECM density colorbar
    cbar = plt.colorbar(shrink=0.65,anchor=(0.4,1))#,label='ECM density')
    cbar.set_label('ECM density')
    cbar.ax.yaxis.set_label_position('left')

    # add contour layer
    # cs = plt.contourf(xx_ecm, yy_ecm, plane_anisotropy, cmap="YlGnBu", levels=levels_ecm)

    # voxel_centres = []


    #### Make square for computing the spheroid area
                
    # for j in cell_df.index:
    #     voxel_side = cell_df_init.loc[0, 'radius'] * 2
    #     voxel_centre_j = voxel_centre(cell_df.loc[j, 'position_x'], cell_df.loc[j, 'position_y'],voxel_side)
    #     pos_x = cell_df.loc[j, 'position_x']
    #     pos_y = cell_df.loc[j, 'position_y']
    #     # print(f'voxel centre={voxel_centre_j}, position={pos_x},{pos_y}')
    #     voxel_centres.append(voxel_centre_j)

    #     color = seaborn.color_palette('muted')[3]

    #     #### Define squares
    #     square = Rectangle((voxel_centre_j[0], voxel_centre_j[1]), width=voxel_side, height=voxel_side, facecolor=color, edgecolor='k' )
    #     ax.add_artist(square)
    
    # # Remove doubles
    # voxel_centres = list(set(voxel_centres))

    # # Compute area
    # voxels_number = len(voxel_centres)
    # spheroid_area = voxels_number * (voxel_side**2)    

    palette = [ '#edf8e9','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#005a32']

    for j in cell_df.index:
        
        #### Cell color depends on number of attached cells        
        attached_cells = int(cell_df.loc[j,'attached_cells'])
        max_attached_cells = 6
        if attached_cells>=max_attached_cells:
            attached_cells = max_attached_cells
        
        #### Assign color to cell depending on neighbours
        color = palette[attached_cells]

        # color = seaborn.color_palette('colorblind')[2]

        #### Define the cell circles
        circ = Circle((cell_df.loc[j, 'position_x'], cell_df.loc[j, 'position_y']), radius=cell_df.loc[j, 'radius'] ,facecolor=color, edgecolor='k')
        ax.add_artist(circ)
        # plt.text(cell_df.loc[j, 'position_x'],  cell_df.loc[j,'position_y'],str(int(cell_df.loc[j,'ID']) ),fontsize=8)


    #### add ECM orientation vectors unscaled by anistorpy ###
    plt.quiver(xx, yy, ECM_x, ECM_y,
    pivot='middle', angles='xy', scale_units='inches', scale=3.0, headwidth=0,headlength=0, headaxislength=0)

    #### Labels and title (will need removed for journal - they will be added manually)
    ax.set_xlabel('x [micron]')
    ax.set_ylabel('y [micron]')
    #fig.colorbar(cs, ax=ax)
    #plt.title('t={t}min'.format(t=t))
    plt.title(f'Ribose={ribose}mM, t={t}min')
    plt.suptitle(f'Prolif rate={prolif}, adh={cell_adh}, rep={cell_rep}, max mot speed={max_mot_speed}', y=0.98)


    #### Carefully place the command to make the plot square AFTER the color bar has been added.
    ax.axis('scaled')
    fig.tight_layout()
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xticks(np.arange(-edge,edge+1, step=40))
    plt.yticks(np.arange(-edge,edge+1, step=40))
    plt.ylim(-edge, edge)
    plt.xlim(-edge, edge)

    #### Colorbar for neighbouring cells
    cmap_2 = mpl.colors.ListedColormap(palette,N=max_attached_cells+1)
    norm = mpl.colors.Normalize(vmin=-0.5,vmax=max_attached_cells+0.5)

    cbar_2 = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap_2),orientation='horizontal', label='Neighbouring cells')
    cbar_2.ax.xaxis.set_label_position('top')
    cbar_2.ax.set_xticklabels(['',0, 1,2,3,4,5, '≥6'])
    
    # fig = plt.gcf()
    # fig.patch.set_alpha(0.0)

    # Plot output

    plt.savefig(save_folder + f'cells_t{int(t)}.png', bbox_inches = "tight")

    plt.close()

    return t

def cell_spheroid_animation(ribose, simulation,seed):

    # Data folder
    data_folder = '../ribose_' + str(ribose) + '/output/output_' + str(simulation)  + '_' + str(seed) + '/'
    
    # Save folder
    save_folder = f'../results/animation/animation_rib{ribose}_{simulation}_{seed}/'   

    images = []
    
    files = os.listdir(data_folder)
    files.sort()

    parameters = []

    for i in range(len(files)):
        if not re.search('_ECM\.mat', files[i]):
            continue
        parameters.append([files[i].split('_')[0], data_folder,save_folder])
        
    t = Parallel(n_jobs=-1)(delayed(cell_spheroid_data)(*p) for p in parameters)

    # t = range(0,5001,10)
    for t_i in t:
        images.append(imageio.v2.imread(save_folder + f'cells_t{int(t_i)}.png'))
    
    writer = imageio.get_writer('../results/animation/cells_animation_rib' + str(ribose) + '_' + str(simulation)  + '_' + str(seed) + '.mp4', fps=10)

    for image in images:
        writer.append_data(image)
    
    writer.close()

    # # make the gif
    # imageio.mimsave(save_folder + 'cells_animation.gif',images)

if __name__ == '__main__':
    
    # Figure resolution
    mpl.rcParams['figure.dpi'] = 300
    
    # Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')
    # fig = plt.gcf()
    # fig.patch.set_alpha(0.0)   
    
    
    # Ribose concentrations
    riboses = [0] #,50,200]
    # ribose = ['50']
    
    simulations = [309,311,312,313,314,315] #range(281,294) 
    seed = 0


    for simulation in simulations:

        # # Initiate figure
        # plt.figure()

        for ribose in riboses:
            return_code = os.system(f'mkdir -p ../results/animation/animation_rib{ribose}_{simulation}_{seed}')
            if return_code != 0:
                print(f'Failed with exit code: {return_code}')
                sys.exit()     

            cell_spheroid_animation(ribose,simulation,seed)





    


