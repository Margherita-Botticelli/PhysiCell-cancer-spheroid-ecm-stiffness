from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os, re
from scipy.interpolate import splrep, splev


def average_total_speeds_data(snapshot,data_folder):
    
    ### Function for average cells speeds over time 
    
    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', data_folder)

    # time
    t = ( mcds.get_time() )

    # Cells
    cell_df = mcds.get_cell_df()
    
    # Get each cell's total speed
    total_speeds = []

    for j in cell_df.index:
        total_speeds.append(cell_df.loc[j,'total_speed'] )

    # Compute average speed
    average_total_speeds = np.mean(total_speeds)
    
    return t, average_total_speeds

def average_total_speeds_plot(data_folder,ribose):
    t = []
    average_total_speeds = []
    
    files = os.listdir(data_folder)
    
    for i in range(len(files)):
        if not re.search('_ECM\.mat', files[i]):
            continue
        
        # Call function to get time and average total speeds values
        t_i, average_total_speeds_i = average_total_speeds_data(files[i].split('_')[0], data_folder)
        
        t.append(t_i)
        average_total_speeds.append(average_total_speeds_i)
          
    # Initiate figure, single plot
    #plt.figure()
    
    # Plot average total speeds over time
    data_plot = plt.plot(t,average_total_speeds,label='{ribose}mM'.format(ribose=ribose))
    
    # Plot smooth function
    color = data_plot[0].get_color()
    bspl = splrep(t,average_total_speeds,k=4,s=5)
    y_smooth = splev(t,bspl)
    plt.plot(t,y_smooth,color=color)
    
    #  Set axis labels
    plt.ylabel("$speed$ [micron/min]",fontsize=15)
    plt.xlabel("t [min]",fontsize=15)
    
    # Set axis ticks
    plt.yticks(np.arange(0,1.1,0.25),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/5),fontsize=13)
    #plt.xticks(np.arange(0, 2000+1,200),fontsize=13)
    
    # Set title, single plot
    #plt.title('Average cells total speeds, ribose {ribose}mM'.format(ribose=ribose))
    
    # Set title, overlaied plots
    plt.title('Average cells total speeds')
    
    plt.legend()
    #plt.show()

def average_mot_vs_mech_speeds_data(snapshot,data_folder):
    
    ### Function for average cells speeds over time 
    
    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', data_folder)

    # time
    t = ( mcds.get_time() )

    # Cells
    cell_df = mcds.get_cell_df()
    
    # Get each cell's speed
    total_speeds = []
    mechanical_speeds = []
    motility_speeds = []

    for j in cell_df.index:
        tot_x = cell_df.loc[j,'total_velocity_x']
        tot_y = cell_df.loc[j,'total_velocity_y']
        tot_z = cell_df.loc[j,'total_velocity_z']
        
        mot_x = cell_df.loc[j,'motility_vector_x']
        mot_y = cell_df.loc[j,'motility_vector_y']
        mot_z = cell_df.loc[j,'motility_vector_z']
        
        tot_speed = np.sqrt((tot_x)**2 + (tot_y)**2 + (tot_z)**2)
        
        mot_speed = np.sqrt((mot_x)**2 + (mot_y)**2 + (mot_z)**2)
        
        mech_speed = np.sqrt((tot_x - mot_x)**2 + (tot_y - mot_y)**2 + (tot_z - mot_z)**2)
        
        if t == 0:
            mot_speed = 0
        
        total_speeds.append(tot_speed)
        mechanical_speeds.append(mech_speed)
        motility_speeds.append(mot_speed)
        #print(cell_df.loc[j,'migration_speed'], motility_speed)
    
    
    # Compute average speeds
    average_total_speeds = np.mean(total_speeds)
    average_mechanical_speeds = np.mean(mechanical_speeds)
    average_motility_speeds = np.mean(motility_speeds)
    #print(average_total_speeds,average_mechanical_speeds,average_motility_speeds)
    
    
    return t, average_total_speeds, average_mechanical_speeds, average_motility_speeds



def average_mot_vs_mech_speeds_plot(data_folder,ribose):
    t = []
    average_total_speeds = []
    average_motility_speeds = []
    average_mechanical_speeds = []
    
    files = os.listdir(data_folder)
    
    for i in range(len(files)):
        if not re.search('_ECM\.mat', files[i]):
            continue
        
        # Call function to get time and average total speeds values
        t_i, average_total_speeds_i, average_mechanical_speeds_i, average_motility_speeds_i = average_mot_vs_mech_speeds_data(files[i].split('_')[0], data_folder)
        
        t.append(t_i)
        
        average_total_speeds.append(average_total_speeds_i)
        average_mechanical_speeds.append(average_mechanical_speeds_i)
        average_motility_speeds.append(average_motility_speeds_i)

    
    # Initiate figure, single plot
    #plt.figure()
    
    # Plot average migration speeds over time
    mot_plot = plt.plot(t,average_motility_speeds,label=f'Mot speed {ribose}mM')
    
    # Plot average migration speeds over time
    color = mot_plot[0].get_color()
    plt.plot(t,average_mechanical_speeds, color=color, alpha=0.5, label=f'Mech speed {ribose}mM')
    
    # # Plot smooth function
    # color = data_plot[0].get_color()
    # bspl = splrep(t,average_total_speeds,k=4,s=5)
    # y_smooth = splev(t,bspl)
    # plt.plot(t,y_smooth,color=color)
    
    #  Set axis labels
    plt.ylabel("speed [micron/min]",fontsize=15)
    plt.xlabel("t [min]",fontsize=15)
    
    # Set axis ticks
    plt.yticks(np.arange(0,1.6,0.25),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,t[-1]/5),fontsize=13)
    #plt.xticks(np.arange(0, 2000+1,200),fontsize=13)
    
    # Set title, single plot
    #plt.title('Average cells total speeds, ribose {ribose}mM'.format(ribose=ribose))
    
    # Set title, overlaied plots
    plt.title('Average cells motility and mechanical speeds')
    
    plt.legend()
    #plt.show()


if __name__ == '__main__':
    
    # Figure resolution
    mpl.rcParams['figure.dpi'] = 300
    
    # Plot style
    plt.style.use('ggplot')
    plt.style.use('seaborn-v0_8-colorblind')
    
    # Initiate figure
    plt.figure()
    
    #ribose = ['0','25','50','100','150','200']
    ribose = ['0','50','200']
    ribose = ['0']
    
    folder = '\\proliferation'
    folder = '\\no_proliferation'
    #folder = '\\output\\'
    
    for i in ribose :
        ribose = i
    
        
        # Data folder
        data_folder = '..\\output\\ribose' + ribose + folder + '\\output\\'
        
        # Save folder
        save_folder = '..\\results\\ribose' + ribose + folder + '\\results\\'  
        
        # Total speeds plot w.r.t. time
        average_total_speeds_plot(data_folder,ribose)
        
        # Total speeds plot w.r.t. time
        #average_mot_vs_mech_speeds_plot(data_folder,ribose)
    
    # fig = plt.gcf()
    # fig.patch.set_alpha(0.0)