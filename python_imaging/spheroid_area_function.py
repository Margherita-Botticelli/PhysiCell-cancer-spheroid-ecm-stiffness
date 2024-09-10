from turtle import color
from pyMCDS_ECM import *
import numpy as np
import os
from skimage import io, draw, util
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar




def spheroid_area_function(data,save_folder='../results/',figure=False):
    
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
    r_anisotropy = data['anisotropy_increase_rate'].iloc[0]
    r_density = data['ecm_density_rate'].iloc[0]
    r_orientation = data['fiber_realignment_rate'].iloc[0]

    # os.system(f'mkdir -p ../images/images_rib{ribose}_{simulation}_{seed}')

    # save_folder = f'../images/images_rib{ribose}_{simulation}_{seed}/' 

    spheroid_area_mask = np.zeros((5000, 5000))
        
    # for each x, y pair
    for pos_x, pos_y, r in zip(position_x, position_y, radius):

        pos_x = (pos_x + 500) * 5
        pos_y = (500 - pos_y) * 5
        # print(pos_x, pos_y, flush=True)
        # the draw library returns a truth array of which points would be the shape of
        # a circle at that point (documentation: https://scikit-image.org/docs/stable/api/skimage.draw.html#skimage.draw.disk)
        # we need to give it the masks shape to generate the correct shape truth array
        circle_coordinates = draw.disk((pos_y, pos_x), r*5, shape=spheroid_area_mask.shape)

        spheroid_area_mask[circle_coordinates] = 1.0

    spheroid_area = np.count_nonzero(spheroid_area_mask)/25

    if figure==True:

        # save, the filename then the object, the util function forces it to be of the format
        # that it's happy tot sa v
        # io.imsave(save_folder + f'mask_{int(t)}.png', util.img_as_ubyte(spheroid_area_mask))


        fig = plt.figure()
        ax = plt.gca()

        cmap = mpl.colors.LinearSegmentedColormap.from_list("", ["black",'white'])
        plt.imshow(spheroid_area_mask,extent=(-500, 500, -500, 500),cmap=cmap)
        

        ax.set_facecolor("black")

        plt.grid(False)

        edge = 500
        plt.xticks([])#(np.arange(-edge,edge+1, step=100),fontsize=12,rotation=45)
        plt.yticks([])#(np.arange(-edge,edge+1, step=100),fontsize=12)
        
        plt.ylim(-edge, edge)
        plt.xlim(-edge, edge)

        scalebar = AnchoredSizeBar(ax.transData,
            100, r'100 [$\mu$m]', 'lower right', 
            pad=0.1,
            color='white',
            frameon=False,
            size_vertical=1,
            fontproperties={'size':15})
        ax.add_artist(scalebar)

        ax.axis('scaled')



        # #### Set title, overlaied plots
        # # plt.suptitle('Clusters and outliers', fontsize = 12)
        # plt.title(f'{prolif=}, {cell_adh=}, {cell_rep=}, {max_mot_speed=}\n{r_density=}, {r_orientation=}, {r_anisotropy=}\nt={int(t)}',fontsize=10)

        #### Plot style
        plt.style.use('ggplot')
        plt.style.use('seaborn-v0_8-colorblind')

        plt.savefig(save_folder + f'/statistics/cell_image_rib{ribose}_{simulation}_{seed}_t{int(t)}.png',bbox_inches='tight')
        
        plt.close()

    return spheroid_area
