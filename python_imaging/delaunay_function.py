from scipy.spatial import Delaunay, distance
from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
from sklearn.cluster import AgglomerativeClustering as AC
from skimage import io, draw, util, color
import seaborn

def delaunay_distance_function(data,save_folder,figure=False):

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
        plt.figure()
        ax = plt.gca()

        plt.grid(False)

        #### Set axes
        ax.set_xlim(-510, 510)
        ax.set_ylim(-510, 510)
        ax.set_xticks(np.arange(-500,501,100))
        ax.set_yticks(np.arange(-500,501,100))

        #### Plot network
        plt.triplot(coordinates[:,0], coordinates[:,1], delaunay_network.simplices)
        plt.plot(coordinates[:,0], coordinates[:,1], 'o')

        #### Plot style
        plt.style.use('ggplot')
        plt.style.use('seaborn-v0_8-colorblind')

        plt.savefig(save_folder + f'clusters/delaunay_rib{ribose}_{simulation}_{seed}_t{int(t)}.png', dpi=600)
        
        plt.close()

    return edges_lengths_mean