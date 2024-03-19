from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, re,math
import imageio


def hist_data(snapshot, data_folder):
    
    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', data_folder)

    # time
    t = ( mcds.get_time() )

    # Cells
    cell_df = mcds.get_cell_df()
    
    dist = []

    for j in cell_df.index:
        # Get data from cell
        position_x = cell_df.loc[j, 'position_x']
        position_y = cell_df.loc[j, 'position_y']
        dist_temp = np.sqrt(position_x**2 + position_y**2)
        
        # Append data to lists
        dist.append(dist_temp)
       
    number_of_cells = len(cell_df.index)
    return t, dist, number_of_cells

def hist_plot(snapshot,data_folder,ribose):
    dist = []
    t=[]
    
    t, dist, number_of_cells = hist_data(snapshot, data_folder)
        
    max_dist = max(dist)
    step = 10
    x_max = int (math.ceil(max_dist / step)) * step
    print(x_max)
    bins = int(x_max / step)
    print(bins)

    # Initiate figure
    plt.figure()

    plt.close()
    plt.hist(dist, range=(0,x_max),bins=bins,rwidth=0.8)

    #Set title
    plt.title(f"Ribose {ribose} mM, Time t={t}min")
    
    # Set ticks
    plt.xticks(np.arange(0,x_max+1,step))
    #plt.yticks(np.arange(0,number_of_cells+1,5),fontsize=13)
    
    #  Set axis labels
    plt.ylabel("Number of cells",fontsize=15)
    plt.xlabel("Distance",fontsize=15)

    
def hist_animation(data_folder,save_folder,ribose):
    dist = []
    t=[]
    
    files = os.listdir(data_folder)
    
    for i in range(len(files)):
        if not re.search('_ECM\.mat', files[i]):
            continue
        t_i, dist_i, number_of_cells = (hist_data(files[i].split('_')[0], data_folder))
        t.append(t_i)
        dist.append(dist_i)
        
    max_dist = max(max(x) for x in dist)

    x_max = math.ceil(max_dist*5)/5

    bins = int(x_max*5)

    # Initiate figure
    plt.figure()
    images = []


    for i in range(0,len(t)):
        plt.close()
        plt.hist(dist[i], range=(0,x_max),bins=bins,rwidth=0.8)
    
        #Set title
        plt.title("Ribose {ribose} mM, Time t {t}min".format(ribose=ribose, t=t[i]))
        
        # Set axis ticks
        plt.yticks(np.arange(0,number_of_cells+1,5),fontsize=13)
        
        #  Set axis labels
        plt.ylabel("Number of cells",fontsize=15)
        plt.xlabel("Distance",fontsize=15)
        
        plt.savefig(save_folder + 'hist_t{t}.png'.format(t=int(t[i])))
        images.append(imageio.imread(save_folder + 'hist_t{t}.png'.format(t=int(t[i]))))
    
    # make the gif
    imageio.mimsave(save_folder + 'hist_animation.gif',images)
    
    
if __name__ == '__main__':
    
    # Figure resolution
    mpl.rcParams['figure.dpi'] = 300
    
    # Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')

    # Ribose concentration    
    ribose = '200'
    
    # Folders
    #folder = '\\proliferation'
    folder = '\\no_proliferation'
    
    # Data folder
    data_folder = '..\\output\\ribose' + ribose + folder + '\\output\\'
    
    # Save folder
    save_folder = '..\\results\\ribose' + ribose + folder +'\\results\\'    
    # Number of cells within radius
    snapshot = 'output00000180'
    hist_plot(snapshot,data_folder,ribose)
    
    
    # Number of cells distribution histogram
    #hist_animation(data_folder,save_folder,ribose)
    

