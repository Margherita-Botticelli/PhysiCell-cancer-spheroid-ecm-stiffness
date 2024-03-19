from scipy import cluster
from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance_matrix
from sklearn.cluster import AgglomerativeClustering as AC
from skimage import io, draw, util, color
import seaborn

def cluster_function(data,save_folder):

    #################### CLUSTERING ALGORITHM #############################

    plt.figure()

    #### Collect the data from the dataframe
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    # seed = data['seed'].iloc[0]
    prolif = data['prolif'].iloc[0]
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    
    seeds = np.unique(data['seed'])
    print(f'{seeds=}',flush=True)

    times = np.unique(data['t'])

    main_cluster_cell_percentages = []

    for seed in seeds:
        for t in times:
            t = int(t)
            # print(f'Time: {t}',flush=True)
            position_x = np.array(data[(data['t'] == t) & (data['seed'] == seed)]['position_x'].iloc[0])
            position_y = np.array(data[(data['t'] == t) & (data['seed'] == seed)]['position_y'].iloc[0])
            radius = np.array(data[(data['t'] == t) & (data['seed'] == seed)]['radius'].iloc[0])

            # print(f'{radius=}',flush=True)

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
            
            #### Percentage of cells in main cluster
            main_cluster_label = unique[np.argmax(unique_counts)]
            main_cluster_cell_count = np.count_nonzero(cluster_labels == main_cluster_label)
            main_cluster_cell_percentage = round(main_cluster_cell_count / n * 100, 2)

            main_cluster_cell_percentages.append(main_cluster_cell_percentage)
            

            if t == 5000:
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

                # io.imsave(save_folder + f'clusters_rib{ribose}_{simulation}_{seed}_t{t}.png', util.img_as_ubyte(cluster_image), check_contrast=False)

                plt.savefig(f'../clusters/clusters_rib{ribose}_{simulation}_{seed}_t{t}.png', dpi=600)
                
                plt.close()
        
    #### Make plot for mean main cluster cell percentage
    main_cluster_cell_percentages = np.reshape(np.array(main_cluster_cell_percentages), (-1,len(times))).astype(float)        
    main_cluster_cell_percentage_mean = np.mean(main_cluster_cell_percentages, axis=0)

    plt.figure()

    plt.plot(times, main_cluster_cell_percentage_mean)

    ### Set axis labels
    plt.ylabel("Percentage",fontsize=15)
    plt.xlabel("Time [min]",fontsize=15)

    ### Set axis ticks
    plt.xticks(np.arange(0, times[-1]+1,times[-1]/5),fontsize=13)

    #### Set title
    plt.suptitle('Percentage of cells in main cluster', fontsize = 12)
    plt.title(f'Prolif={prolif}, adh={cell_adh}, rep={cell_rep}, {max_mot_speed=}, t={int(t)}',fontsize=10)

    #### Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')

    plt.savefig(f'../results/main_cluster_cell_percentage_rib{ribose}_{simulation}.png', dpi=600)
            
    plt.close()
    



        












    # #### Array with each cell's coordinates
    # X = np.vstack((position_x[0][time],position_y[0][time])).T
    # # print(X,flush=True)

    # #### Set minimum distance between points 
    # radius = radii[0][0]
    # eps = 1.25 * radius * 2
    # # eps=20

    # #### Set minimum number of elements per cluster
    # min_samples = 5

    # # #### Find clusters
    # # dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    # # dbscan.fit(X)

    # # #### Labels with cluster number (if -1 they are outliers)
    # # labels = dbscan.labels_
    # # # print(labels,flush=True)

    # # #### Core points indices
    # core_points_indices = []
    # # core_points_indices = dbscan.core_sample_indices_
    # # print(core_points_indices,flush=True)


    # # Initialize the OPTICS clustering model
    # optics = OPTICS(min_samples=min_samples, xi=0.05, min_cluster_size=0.05)
    
    # # Fit the model to the data
    # optics.fit(X)
    
    # # Get the cluster labels
    # labels = optics.labels_
    
    # #### Identify the outliers
    # outliers = np.where(labels == -1)[0]

    # #### Set size of each dot
    # r = [x for x in radii[time]]
    # r_outliers = [r[o] for o in outliers]

    # #### Color list
    # color_list = seaborn.color_palette('colorblind') + seaborn.color_palette('dark') + seaborn.color_palette('muted') + seaborn.color_palette('bright') + ['black']
    # colors = [color_list[l] for l in labels]

    # for i in range(len(r)):
    #     if(i in core_points_indices):
    #         edgecolor = 'white'
    #     else:
    #         edgecolor = 'black'
            
    #     #define circles
    #     circle = plt.Circle((X[i, 0], X[i, 1]), radius=r[i], facecolor=colors[i], edgecolor=edgecolor)

    #     #add circles to plot
    #     plt.gca().add_artist(circle)


    # # #### Plot the data with the outliers highlighted
    # # plt.scatter(X[:, 0], X[:, 1], c=colors, edgecolors='black', s=s)
    # # plt.scatter(X[outliers, 0], X[outliers, 1], c="black", s=s_outliers)

    # ### Set axis ticks
    # plt.xticks(range(-500,500,40), fontsize=10)
    # plt.yticks(range(-500,500,40),fontsize=10)
    # plt.xlim(-400,400)
    # plt.ylim(-400,400)
    # plt.xticks(rotation=45, ha="right")

    # ### Make axis same size
    # plt.axis('scaled')

    
    # #### Set title, overlaied plots
    # plt.title('Clusters and outliers')
    # plt.suptitle(f'Prolif={prolif}, adh={cell_adh}, rep={cell_rep}, max mot speed={max_mot_speed}, t={t_current}', y=0.98)
    
    # #### Save figure
    # plt.savefig(save_folder + 'cluster_outliers_rib' + str(ribose)+ '_' + str(simulation) + '.png')#, bbox_inches = "tight")
