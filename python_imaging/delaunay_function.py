from scipy.spatial import Delaunay, distance
from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
from sklearn.cluster import AgglomerativeClustering as AC
from skimage import io, draw, util, color
import seaborn
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.patches import Rectangle

def delaunay_distance_function(data,save_folder='../results/',figure=False):

    #################### CLUSTERING ALGORITHM #############################

    #### Collect the data from the dataframe
    ribose = data['ribose'].iloc[0]
    simulation = data['simulation'].iloc[0]
    seed = data['seed'].iloc[0]
    t = data['t'].iloc[0]
    radius = data['radius'].to_numpy()
    position_x = data['position_x'].to_numpy()
    position_y = data['position_y'].to_numpy()
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    # r_anisotropy = data['anisotropy_increase_rate'].iloc[0]
    r_density = data['ecm_density_rate'].iloc[0]
    # r_orientation = data['fiber_realignment_rate'].iloc[0]

    #### Number of cells
    n = len(position_x)
        
    #### Find distance between cells centres
    coordinates = np.stack((position_x,position_y), axis=1)
    delaunay_network = Delaunay(coordinates) 

    #### Find edges 
    edges = []

    #### Find edges in all simplices
    for tri in delaunay_network.simplices:
        #### Append point indices
        edges.append([tri[0],tri[1]])
        edges.append([tri[1],tri[2]])
        edges.append([tri[0],tri[2]])
    edges = np.array(edges)
    edges = np.sort(edges, axis=1)
    # np.set_printoptions(threshold=sys.maxsize)
    # print(f'{edges=}, {len(edges)}',flush=True)

    #### Remove duplicate edges (sorting necessary for unique to work)
    edges_unique = np.unique(edges,axis=0)
    # print(f'{edges_unique=}, {len(edges_unique)},',flush=True)

    #### Find point coordinates for the edges
    edges_unique_coords = coordinates[edges_unique].astype(np.float64)
    # print(f'{np.min(edges_unique_coords)=}, {np.max(edges_unique_coords)=}',flush=True)

    # print(f'{edges_unique_coords=}',flush=True)

    #### Find edges lengths
    edges_lengths = []
    for p1, p2 in edges_unique_coords:
        dist = np.linalg.norm(p1-p2)
        if dist == np.inf:
            print(f'{p1=}, {p2=}, {dist=}', flush=True)
        edges_lengths.append(dist)

    #### Compute means
    edges_lengths_mean = np.mean(edges_lengths)
    # print(f'{simulation=}, {edges_lengths_mean=}',flush=True)

    if figure == True:
        #### Plot network
        #### Initialise figure
        plt.figure(figsize=(5,5))
        ax = plt.gca()

        plt.grid(False)

        #### Set axes

        edge = 500
        plt.xticks([])#(np.arange(-edge,edge+1, step=100),fontsize=12,rotation=45)
        plt.yticks([])#(np.arange(-edge,edge+1, step=100),fontsize=12)

        plt.ylim(-edge, edge)
        plt.xlim(-edge, edge)
    
        #### Plot network
        plt.plot(coordinates[:,0], coordinates[:,1], 'o',color=seaborn.color_palette('colorblind')[2],zorder=0)
        plt.triplot(coordinates[:,0], coordinates[:,1], delaunay_network.simplices, color='black',zorder=1,lw=0.5)

        # ax.add_patch( Rectangle((-500, -500), 1000, 1000, fc='none', ec='black') )

        scalebar = AnchoredSizeBar(ax.transData,
                    100, r'100 [$\mu$m]', 'lower right', 
                    pad=0.1,
                    color='black',
                    frameon=False,
                    size_vertical=1,
                    fontproperties={'size':15})
        ax.add_artist(scalebar)




        # #### Plot style
        # plt.style.use('ggplot')
        # plt.style.use('seaborn-v0_8-colorblind')

        plt.savefig(save_folder + f'/statistics/delaunay_rib{ribose}_{simulation}_{seed}_t{int(t)}.png',bbox_inches='tight')
        
        plt.close()

    return edges_lengths_mean