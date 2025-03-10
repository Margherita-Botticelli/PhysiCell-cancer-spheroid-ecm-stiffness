from scipy.spatial import Delaunay
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.patches import Rectangle
import seaborn

def delaunay_distance_function(data, save_folder='../results/', figure=False):
    """
    Calculate the mean Delaunay distance and optionally plot the Delaunay network.

    Parameters:
    - data: pandas DataFrame containing the simulation data.
    - save_folder: Directory path to save the plot (default is '../results/').
    - figure: Boolean to decide whether to generate and save a figure (default is False).

    Returns:
    - edges_lengths_mean: Mean length of the edges in the Delaunay network.
    """
    
    #### Extract simulation parameters from the DataFrame
    # ribose = data['ribose'].iloc[0]
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
    r_density = data['ecm_density_rate'].iloc[0]

    #### Number of cells
    n = len(position_x)
    
    #### Prepare coordinates for Delaunay triangulation
    coordinates = np.stack((position_x, position_y), axis=1)
    delaunay_network = Delaunay(coordinates)

    #### Extract edges from the Delaunay triangulation
    edges = []
    for tri in delaunay_network.simplices:
        edges.append([tri[0], tri[1]])
        edges.append([tri[1], tri[2]])
        edges.append([tri[0], tri[2]])
    edges = np.array(edges)
    edges = np.sort(edges, axis=1)

    #### Remove duplicate edges
    edges_unique = np.unique(edges, axis=0)

    #### Calculate edge lengths
    edges_unique_coords = coordinates[edges_unique].astype(np.float64)
    edges_lengths = [np.linalg.norm(p1 - p2) for p1, p2 in edges_unique_coords]

    #### Compute mean edge length
    edges_lengths_mean = np.mean(edges_lengths)

    #### Optionally plot the Delaunay network
    if figure:
        plt.figure(figsize=(5, 5))
        ax = plt.gca()

        plt.grid(False)
        edge = 500
        plt.xticks([])  # Remove x-axis ticks
        plt.yticks([])  # Remove y-axis ticks
        plt.ylim(-edge, edge)
        plt.xlim(-edge, edge)

        #### Plot Delaunay network
        plt.plot(coordinates[:, 0], coordinates[:, 1], 'o', color=seaborn.color_palette('colorblind')[2], zorder=0)
        plt.triplot(coordinates[:, 0], coordinates[:, 1], delaunay_network.simplices, color='black', zorder=1, lw=0.5)

        #### Add a scale bar to the plot
        scalebar = AnchoredSizeBar(ax.transData, 100, r'100 [$\mu$m]', 'lower right', 
                    pad=0.1, color='black', frameon=False, size_vertical=1, 
                    fontproperties={'size': 15})
        ax.add_artist(scalebar)

        #### Save the plot
        # plt.savefig(save_folder + f'/statistics/delaunay_rib{ribose}_{simulation}_{seed}_t{int(t)}.png', bbox_inches='tight')
        plt.savefig(save_folder + f'/statistics/delaunay_{simulation}_{seed}_t{int(t)}.png', bbox_inches='tight')
        plt.close()

    return edges_lengths_mean
