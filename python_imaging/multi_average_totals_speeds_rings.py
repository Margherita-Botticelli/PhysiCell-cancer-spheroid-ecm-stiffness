from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, re
from scipy.interpolate import splrep, splev

    
def average_total_speeds_rings_data(snapshot,data_folder,R):
    
    ### Function for average cells speeds over time 
    
    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', data_folder)

    # time
    t = ( mcds.get_time() )

    # Cells
    cell_df = mcds.get_cell_df()
    
    # Get each cell's position and total speed
    positions_x = []
    positions_y = []
    total_speeds = []
    dist = []
    
    for j in cell_df.index:
        # Get data from cell
        position_x = cell_df.loc[j, 'position_x']
        position_y = cell_df.loc[j, 'position_y']
        dist_temp = np.sqrt(position_x**2 + position_y**2)
        
        # Append data to lists
        positions_x.append(position_x)
        positions_y.append(position_y)
        total_speeds.append(cell_df.loc[j,'total_speed'])
        dist.append(dist_temp)
       
    # Turn lists into arrays
    dist = np.array(dist)
    total_speeds = np.array(total_speeds)
    
    
    # Split into rings
    total_speeds_R = [[] for _ in range(len(R))]

    for i in range(0, len(dist)):
        over = True
        for j in range(0,len(R)-1):
            if R[j] <= dist[i] < R[j+1]:
                total_speeds_R[j].append(total_speeds[i])
                over = False
        if over == True:
            total_speeds_R[len(R)-1].append(total_speeds[i])
    
    # Compute average speed
    average_total_speeds_R =  [[] for _ in range(len(R))]
    
    for j in range(0,len(R)):
        average_total_speeds_R[j] = np.mean(total_speeds_R[j])

    
    #print(total_speeds_Rover)
    
    return t, average_total_speeds_R



def average_total_speeds_rings_plot(data_folder,ribose,R):

    t = []
    average_total_speeds_R = [[] for _ in range(len(R))]
    
    files = os.listdir(data_folder)
    
    for i in range(len(files)):
        if not re.search('_ECM\.mat', files[i]):
            continue
        
        # Call function to get time and average total speeds values
        t_i, average_total_speeds_R_i = average_total_speeds_rings_data(files[i].split('_')[0], data_folder,R)
    
        
        t.append(t_i)
        for j in range(0,len(R)):
            average_total_speeds_R[j].append(average_total_speeds_R_i[j])

    
    # Initiate figure
    plt.figure()
    
    #R_plot =  [[] for _ in range(len(R))]
    
    # Plot average total speeds over time
    for j in range(0,len(R)-1):
        plt.plot(t,average_total_speeds_R[j],label=f"[R{j},R{j+1})")
    plt.plot(t,average_total_speeds_R[len(R)-1],label=r"$\geq R{n}$".format(n=len(R)-1))

    # Legend
    plt.legend()
    
    # # Plot smooth function
    # color = R_plot[j][0].get_color()
    # bspl = splrep(t, average_total_speeds_R0, k=3, s=5)
    # y_smooth = splev(t, bspl)
    # plt.plot(t, y_smooth, color=color)
    
    #  Set axis labels
    plt.ylabel("speed [micron/min]",fontsize=15)
    plt.xlabel("t [min]",fontsize=15)
    
    # Set axis ticks
    plt.yticks(np.arange(0,1.6,0.25),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/5),fontsize=13)
    
    r_display = [f' R{n}={v}' for n, v in enumerate(R)]
    title = f'Average cells total speeds, ribose {ribose}mM\nRadii'
    
    for s in r_display:
        title = title + s 
    
    # Set title
    plt.title(title)
     


if __name__ == '__main__':
    
    # Figure resolution
    mpl.rcParams['figure.dpi'] = 300
    
    # Plot style
    #plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')
    

    ribose = '200'
    
    folder = '\\proliferation'
    folder = '\\no_proliferation'
    #folder = '\\output\\'
    
    # Data folder
    data_folder = '..\\output\\ribose' + ribose + folder + '\\output\\'
    
    # Save folder
    save_folder = '..\\results\\ribose' + ribose + folder + '\\results\\'  
    # Make rings list
    step = 15
    R_min = 0
    R_max = 60
    
    R = np.arange(R_min, R_max+1,step)
    
    # Total speeds plot w.r.t. time
    average_total_speeds_rings_plot(data_folder,ribose,R)

    # fig = plt.gcf()
    # fig.patch.set_alpha(0.0)
    
    
