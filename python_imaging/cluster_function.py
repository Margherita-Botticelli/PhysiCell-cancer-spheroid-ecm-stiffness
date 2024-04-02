from scipy import cluster
from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
from sklearn.cluster import AgglomerativeClustering as AC
from skimage import io, draw, util, color
import seaborn

def cluster_function(data,save_folder,figure=False):

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
    distances = distance_matrix(coordinates,coordinates) 

    #### Compute sum of radii
    radii = np.repeat(radius,n).reshape(n,n)
    radii = radii + radii.T
    # print(f'radii:{radii}',flush=True)
    radii = radii * 1.25

    #### Find distance between cells surface 
    interaction_distances = distances - radii   # If the interaction distance is positive the cells are not interacting

    #### If the interaction distance is greater than zero, then the two cells are not interacting
    cluster_metric = AC(metric='precomputed', linkage = 'single', distance_threshold= 0.0, n_clusters=None)

    clusters = cluster_metric.fit(interaction_distances)

    cluster_labels = clusters.labels_

    #### Find outliers and label them as -2
    unique, unique_counts = np.unique(cluster_labels,return_counts=True)
    # print(f'{unique=} and {unique_counts=}',flush=True)
    
    for u,count in zip(unique,unique_counts):
        if count < 3:
            cluster_labels[cluster_labels == u] = -2

    #### Clustering metrics
    #### Number of clusters ignoring outliers
    number_clusters = len(np.unique(cluster_labels)) -1
    # print(f'{number_clusters=}',flush=True)
    
    
    #### Percentage of cells in main cluster
    main_cluster_label = unique[np.argmax(unique_counts)]
    main_cluster_cell_count = np.count_nonzero(cluster_labels == main_cluster_label)
    main_cluster_cell_percentage = round(main_cluster_cell_count / n * 100, 2)
    # print(f'{main_cluster_cell_percentage=}',flush=True)


    if figure == True and seed==0:
        ### Make plot with different colours
        cluster_mask = np.zeros((5000, 5000))
            
        # for each x, y pair
        for pos_x, pos_y, r, label in zip(position_x, position_y, radius, cluster_labels):
            pos_x = (pos_x + 500) * 5
            pos_y = (500 - pos_y) * 5
            # print(pos_x, pos_y, flush=True)
            circle_coordinates = draw.disk((pos_y, pos_x), r*5, shape=cluster_mask.shape)

            cluster_mask[circle_coordinates] = label +1

        # color_list = [(0.0, 0.0, 0.0)]
        color_list = []
        color_list.extend(list(seaborn.color_palette('colorblind')))
        color_list.extend(list(seaborn.color_palette('dark')))
        color_list.extend(list(seaborn.color_palette('muted')))
        color_list.extend(list(seaborn.color_palette('bright')))

        cluster_image = color.label2rgb(cluster_mask,bg_color='white', colors=color_list)
        cluster_image[cluster_mask == -1] = (0.0,0.0,0.0)
        
        plt.figure()
        ax = plt.gca()

        plt.imshow(cluster_image,extent=(-500, 500, -500, 500))

        plt.grid(False)

        ax.set_xlim(-510, 510)
        ax.set_ylim(-510, 510)
        ax.set_xticks(np.arange(-500,501,100))
        ax.set_yticks(np.arange(-500,501,100))

        #### Set title, overlaied plots
        # plt.suptitle('Clusters and outliers', fontsize = 12)
        plt.title(f'Prolif={prolif}, adh={cell_adh}, rep={cell_rep}, {max_mot_speed=}, t={int(t)}\n{number_clusters=}, {main_cluster_cell_percentage=}',fontsize=10)

        #### Plot style
        plt.style.use('ggplot')
        plt.style.use('seaborn-v0_8-colorblind')

        plt.savefig(save_folder + f'clusters/clusters_rib{ribose}_{simulation}_{seed}_t{int(t)}.png', dpi=600)
        
        plt.close()

    return main_cluster_cell_percentage