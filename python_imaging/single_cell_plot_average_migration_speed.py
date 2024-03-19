from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys



   
def cell_migration_speed_over_time(folder,color):
    
    ### Function for plot of cell speed over time
    
    # Open migration speed file, taken at each mechanics time step
    with open(folder + 'migration_speed.txt') as f:
        ms = f.read()
    ms = ms.split('\n')
    ms = [float(n) for n in ms if n]
    print(ms)
    
    # Load cell and microenvironment data
    mcds = pyMCDS('initial.xml', folder)  
    
    # Cells
    cell_df = mcds.get_cell_df()
    ribose_concentration = float(cell_df['ribose_concentration'])
    
    # Compute average speed
    average_speed = np.mean(ms)
    
    # Time spacing
    t = np.arange(0,len(ms)*0.1,0.1)
    
    ###### Plot migration speed for every time step ######
    
    # Initiate figure
    plt.figure()
    fig1, ax = plt.subplots()
    
    # Set axis and box aspect ration
    ax.set_xlim(0 , 1)
    ax.set_box_aspect(0.7)

    # Plot of speed over time
    plt.plot(t,ms,color=color,label='ribose ' + r'$\zeta_0={n}mM$'.format(n=int(ribose_concentration)))
    
    # Horizontal line for average speed
    plt.hlines(y=average_speed, xmin=0,xmax=t[-1],color=list(plt.rcParams['axes.prop_cycle'])[3]['color'])
    
    # Box with average speed value
    bbox = dict(boxstyle='round',edgecolor=list(plt.rcParams['axes.prop_cycle'])[3]['color'],facecolor = 'white', alpha = 0.7)
    plt.text(t[-1]/2,average_speed, f'{average_speed:.3f}',color=list(plt.rcParams['axes.prop_cycle'])[3]['color'], horizontalalignment='center',fontsize=15, bbox=bbox)
    
    #  Set axis labels
    plt.ylabel("$s_{mot}$ [micron/min]",fontsize=15)
    plt.xlabel("t [min]",fontsize=15)
    
    # Set axis ticks
    plt.yticks(np.arange(0,1.1,0.25),fontsize=13)
    plt.xticks(np.arange(0, t[-1]+1,step=200),fontsize=13)

    # Set legend
    plt.legend(loc='lower center', fontsize=13,facecolor = 'white')
    
    # Plot tifle
    if (ribose_concentration==0):
        plt.title('Cancer cells motility speed',fontsize=18)
    
    # Save plot figure
    #plt.savefig(save_folder + 'migration_speed_plot_average' + '.jpeg')
    #plt.axis('scaled')
    
    # Return average speed
    return average_speed

    



#sys.exit()

if __name__ == '__main__':
    
    # Figure resolution
    mpl.rcParams['figure.dpi'] = 300
    
    # Plot style
    plt.style.use('ggplot')
    
    # Data folder
    folder = '..\\output\\'
    
    # Save folder
    save_folder = '..\\user_projects\\ecm_extension\\results\\'
    
    # Call function for plots of cell speed over time 
    cell_migration_speed_over_time(folder, '#366AB3')
    #cell_speed_over_time('..\\results\\poster\\ribose50\\', '#147A00')
    #cell_speed_over_time('..\\results\\poster\\ribose200\\', '#B86B56')
    
    
    ###### Histogram plot of average speeds with varying ribose concentration ######
    
    # Compute average speed for different ribose concentrations 
    
    # Initialise average_speed array
    average_speed = []
    
    # Folder with data for different ribose concentrations
    folder = ['..\\results\\poster\\ribose0\\','..\\results\\poster\\ribose50\\','..\\results\\poster\\ribose200\\']
    
    # Read migration speed file for each ribose concentration
    for folder in folder:
        with open(folder + 'migration_speed.txt') as f:
            ms = f.read()
        ms = ms.split('\n')
        ms = [float(n) for n in ms if n]
        print(ms)
    
        # Append average speed 
        average_speed.append(np.mean(ms))
 
    # Initiate figure
    plt.figure()
    fig1, ax = plt.subplots()
    
    # Set axis and box aspect ration
    ax.set_ylim(0 ,1.05)
    ax.set_box_aspect(1.5)    
    
    # One bar for each ribose concentration
    hbar0 = plt.bar('0mM', round(average_speed[0],3), color='#366AB3',width=0.8)
    hbar50 = plt.bar('50mM', round(average_speed[1],3), color='#147A00',width=0.8)
    hbar200 = plt.bar('200mM', round(average_speed[2],3), color='#B86B56',width=0.8)
    
    # Labels on top of bars
    plt.bar_label(hbar0,color='#366AB3',fontsize=13)
    plt.bar_label(hbar50,color='#147A00',fontsize=13)
    plt.bar_label(hbar200,color='#B86B56',fontsize=13)
    
    # Set axis ticks
    plt.yticks(np.arange(0,1.1,0.25),fontsize=13)
    plt.xticks(fontsize=13)
    
    # Set axis labels
    plt.xlabel('Ribose ' + r'$\zeta_0$', fontsize=15)#
    plt.ylabel("Average $s_{mot}$",fontsize=15)
    
    # Plot title    
    plt.title('Average motility speed',fontsize=18)
    
    plt.show




# images = []

# for n in range(0,len(t)):

    
#     #plt.figure()
#     plt.plot(t,migration_speed)
#     #plt.plot(t[n],migration_speed[n],'.')
#     plt.vlines(x=t[n], ymin=0,ymax=1,color='orange')
 
#     plt.ylabel("$s_{mot}$ [micron/min]",fontsize=12)
#     plt.yticks(np.arange(0,1.1,0.25),fontsize=10)
 
#     plt.xlabel("t [min]",fontsize=12)
#     step = t[-1] / 10
#     plt.xticks(np.arange(0, t[-1]+1, step),fontsize=10)
#     # ax1.grid(color='r', linestyle='-')
#     plt.title('Cancer cells migration speed\n Ribose concentration={r}mM'.format(r=ribose_concentration),fontsize=15)
    
#     plt.savefig(folder + 'migration_speed_plot_average_{n}'.format(n=n) + '.jpeg')
    
#     images.append(imageio.imread(folder + 'migration_speed_plot_average_{n}'.format(n=n) + '.jpeg'))
#     plt.close()
        
    
# # make the gif
# imageio.mimsave(folder + 'migration_speed_plot_average' + '.gif',images)










