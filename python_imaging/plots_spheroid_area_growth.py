from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn


def plots_spheroid_area_growth(data,fig_num,save_folder,simulation_name):

    #################### PLOT SPHEROID AREA SLOPES #############################

    #### Initiate figure
    plt.figure(fig_num)

    ribose = data['ribose'].iloc[0]
    t = data['t'].iloc[0]
    prolif = data['prolif'].iloc[0]
    cell_adh = float(data['cell_adh'].iloc[0])
    cell_rep = float(data['cell_rep'].iloc[0])
    max_mot_speed = data['max_mot_speed'].iloc[0]
    spheroid_area = data['spheroid_area'].to_numpy()

    spheroid_area = np.vstack(spheroid_area)

    #### Select color for plot depending on ribose concentration
    if(ribose == 0):
        color_rib = seaborn.color_palette('colorblind')[0]
    elif(ribose == 50):
        color_rib = seaborn.color_palette('colorblind')[1]
    elif(ribose == 200):
        color_rib = seaborn.color_palette('colorblind')[2]

    spheroid_area = np.mean(spheroid_area, axis=0)
    
    spheroid_area_init = spheroid_area[0]
    spheroid_area_fin = spheroid_area[-1]

    spheroid_area_ratio = spheroid_area_fin/spheroid_area_init

    # t_slope = t[10:len(t)]

    # res = stats.linregress(t_slope,spheroid_area[10:len(t)])

    # plt.errorbar(f'{max_mot_speed}', spheroid_area_ratio, color=color_rib, marker="o", linestyle="none") #, label=f'{int(label_rib)} mM')
    plt.errorbar(f'({cell_adh},{cell_rep})', spheroid_area_ratio, color=color_rib, marker="o", linestyle="none") #, label=f'{int(label_rib)} mM')

    ###  Set axis labels
    plt.ylabel("Growth",fontsize=15)
    # plt.xlabel("Max motility speed",fontsize=15)
    plt.xlabel("(Adhesion, repulsion)",fontsize=15)
    plt.xticks(rotation=45, ha="right")

    ### Set title, overlaied plots
    plt.title(f'Spheroid growth relative to t$_0$', y=1.0)
    # plt.suptitle(f'Prolif rate:{prolif}, adhesion:{cell_adh}, repulsion:{cell_rep}', y=0.98)
    plt.suptitle(f'Prolif rate:{prolif}, max mot speed:{max_mot_speed}, ratio adhesion/repulsion:{cell_adh/cell_rep}', y=0.98)

    legend_drawn_flag = True
    
    plt.legend(["0 mM", "50 mM", '200 mM'],title='Ribose',frameon=legend_drawn_flag)
    
    # plt.savefig(save_folder + f'plots/spheroid_growth_rib{ribose}_{simulation_name}.png', bbox_inches = "tight")
    plt.savefig(save_folder + f'plots/spheroid_growth_rib{ribose}_{simulation_name}.png', bbox_inches = "tight")
    # return res.slope

