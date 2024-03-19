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

    number_of_cells = len(cell_df.index)
    
    return t, number_of_cells

def hist_plot(snapshot,folder,ribose):

    # Initiate distance array
    n = len(ribose)
    number_of_cells =  [[] for _ in range(n)]
    
    for i in np.arange(0,n):
        
        # Data folder
        data_folder = '..\\output\\ribose' + str(ribose[i]) + folder + '\\output\\'
        
        t, number_of_cells[i] = hist_data(snapshot, data_folder)
    
    print(number_of_cells)
    
    # One bar for each ribose concentration
    hbar0 = plt.bar('0mM', number_of_cells[0],width=0.8)
    hbar50 = plt.bar('50mM', number_of_cells[1], width=0.8)
    hbar200 = plt.bar('200mM', number_of_cells[2], width=0.8)
    
    # # Labels on top of bars
    # plt.bar_label(hbar0,fontsize=15)
    # plt.bar_label(hbar50,fontsize=15)
    # plt.bar_label(hbar200,fontsize=15)
    

    #Set title
    plt.title(f"Total number of cells at time t={t}min")
    
    # Set ticks
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    
    #  Set axis labels
    plt.ylabel("Number of cells",fontsize=15)
    plt.xlabel("Ribose",fontsize=15)

    
# def hist_animation(data_folder,save_folder,ribose):
#     dist = []
#     t=[]
    
#     files = os.listdir(data_folder)
    
#     for i in range(len(files)):
#         if not re.search('_ECM\.mat', files[i]):
#             continue
#         t_i, dist_i, number_of_cells = (hist_data(files[i].split('_')[0], data_folder))
#         t.append(t_i)
#         dist.append(dist_i)
        
#     max_dist = max(max(x) for x in dist)

#     x_max = math.ceil(max_dist*5)/5

#     bins = int(x_max*5)

#     # Initiate figure
#     plt.figure()
#     images = []


#     for i in range(0,len(t)):
#         plt.close()
#         plt.hist(dist[i], range=(0,x_max),bins=bins,rwidth=0.8)
    
#         #Set title
#         plt.title("Ribose {ribose} mM, Time t {t}min".format(ribose=ribose, t=t[i]))
        
#         # Set axis ticks
#         plt.yticks(np.arange(0,number_of_cells+1,5),fontsize=13)
        
#         #  Set axis labels
#         plt.ylabel("Number of cells",fontsize=15)
#         plt.xlabel("Distance",fontsize=15)
        
#         plt.savefig(save_folder + 'hist_t{t}.png'.format(t=int(t[i])))
#         images.append(imageio.imread(save_folder + 'hist_t{t}.png'.format(t=int(t[i]))))

if __name__ == '__main__':
    
    # Figure resolution
    mpl.rcParams['figure.dpi'] = 300
    
    # Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')

    # Ribose concentration    
    ribose = [0,50,200]
    # ribose = [50]
    
    # Folders
    # folder = '\\proliferation'
    folder = '\\proliferation\\mot_speed_1'
    # folder = '\\no_proliferation'
    
    # Snapshot
    snapshot = 'output00000200'
    
    # Call histogram plot function
    hist_plot(snapshot,folder,ribose)
    
    
    # Number of cells distribution histogram
    #hist_animation(data_folder,save_folder,ribose)
    

