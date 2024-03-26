from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn

def cell_velocities(data,save_folder):
    #### Collect the data from the dataframe
    simulation = data['simulation'].iloc[0]
    ribose = data['ribose'].iloc[0]
    # seed = data['seed'].iloc[0]
    prolif = round(float(data['prolif'].iloc[0]),5)
    cell_adh = data['cell_adh'].iloc[0]
    cell_rep = data['cell_rep'].iloc[0]
    max_mot_speed = data['max_mot_speed'].iloc[0]
    
    times = np.unique(data['t'])

    for t in times:
        t = int(t)
        # print(f'Time: {t}',flush=True)
        total_velocity_x = np.array(data[(data['t'] == t)]['total_velocity_x'])
        total_velocity_y = np.array(data[(data['t'] == t)]['total_velocity_y'])
    
        motility_velocity_x = np.array(data[(data['t'] == t)]['motility_vector_x'])
        motility_velocity_y = np.array(data[(data['t'] == t)]['motility_vector_y'])
        # print(f"{motility_velocity_x=}")
        # print(f"{motility_velocity_y=}")
                
        mechanics_velocity_x = total_velocity_x - motility_velocity_x
        mechanics_velocity_y = total_velocity_y - motility_velocity_y
        
        root = np.sqrt((np.square(total_velocity_x) + np.square(total_velocity_y)))
        root[root == 0] = 1
        root = 1
        total_velocity_x_norm = total_velocity_x / root
        total_velocity_y_norm = total_velocity_y / root
        # print(total_velocity_x_norm)

        root = np.sqrt((np.square(motility_velocity_x) + np.square(motility_velocity_y)))
        root[root == 0] = 1
        root = 1
        motility_velocity_x_norm = motility_velocity_x / root
        motility_velocity_y_norm = motility_velocity_y / root
        # print(f"{motility_velocity_y_norm=}")

        root = np.sqrt((np.square(mechanics_velocity_x) + np.square(mechanics_velocity_y)))
        root[root == 0] = 1
        root = 1
        mechanics_velocity_x_norm = mechanics_velocity_x / root
        mechanics_velocity_y_norm = mechanics_velocity_y / root 
        # total_velocity_x_norm[np.isnan(total_velocity_x_norm)] = 0
        # total_velocity_y_norm[np.isnan(total_velocity_y_norm)] = 0
        # motility_velocity_x_norm[np.isnan(motility_velocity_x_norm)] = 0
        # motility_velocity_y_norm[np.isnan(motility_velocity_y_norm)] = 0
        # mechanics_velocity_x_norm[np.isnan(mechanics_velocity_x_norm)] = 0
        # mechanics_velocity_y_norm[np.isnan(mechanics_velocity_y_norm)] = 0

        # print(f"{motility_velocity_y_norm=}")

        fig, axs = plt.subplots(2, 3,sharex=True, sharey=True)
        # #### Plot style
        # fig.style.use('ggplot')
        # fig.style.use('seaborn-v0_8-colorblind')

        fig.suptitle(f'Prolif={prolif}, adh={cell_adh}, rep={cell_rep}, {max_mot_speed=}, t={int(t)}',fontsize=10)
   
        bins =  [-1.1,-0.9,-0.7,-0.5,-0.3,-0.1,0.1,.3,.5,.7,.9,1.1]
        axs[0, 0].hist(total_velocity_x_norm,bins=bins)
        #### Set title
        axs[0, 0].set_title('Total velocities x', fontsize = 12)
        axs[0,0].set_ylim(0,100)
        axs[0,0].set_xticks([-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1])

        axs[1, 0].hist(total_velocity_y_norm,bins=bins)
        #### Set title
        axs[1, 0].set_title('Total velocities y', fontsize = 12)
   

        axs[0, 1].hist(motility_velocity_x_norm,bins=bins)
        #### Set title
        axs[0, 1].set_title('Motility velocities x', fontsize = 12)

        axs[1, 1].hist(motility_velocity_y_norm,bins=bins)
        #### Set title
        axs[1, 1].set_title('Motility velocities y', fontsize = 12)

        axs[0, 2].hist(mechanics_velocity_x_norm,bins=bins)
        #### Set title
        axs[0, 2].set_title('Mechanics velocities x', fontsize = 12)

        axs[1, 2].hist(mechanics_velocity_y_norm,bins=bins)
        #### Set title
        axs[1, 2].set_title('Mechanics velocities y', fontsize = 12)

        plt.savefig(save_folder + f'plots/cell_velocities_rib{ribose}_{simulation}_t{t}.png', dpi=600)
        plt.close()




