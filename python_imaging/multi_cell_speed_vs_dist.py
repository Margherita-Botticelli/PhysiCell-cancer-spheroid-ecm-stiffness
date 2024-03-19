from pyMCDS_ECM import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import splrep, splev


def cell_speed_heatmap(snapshot, data_folder, save_folder):
    
    #################################  Load data  ########################
    # load cell and microenvironment data
    mcds = pyMCDS(snapshot + '.xml', data_folder)

    # loads and reads ECM data
    mcds.load_ecm(snapshot + '_ECM.mat', data_folder)

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
    
    # Sort arrays by ordering the dist array 
    indices = np.argsort(dist)
    x = dist[indices]
    y = total_speeds[indices]
    
    
    def movingaverage(arr, window_size):
        moving_averages = []
        i = 0
        while i < len(arr) - window_size + 1:
            
            # Store elements from i to i+window_size
            # in list to get the current window
            window = arr[i : i + window_size]

          
            # Calculate the average of current window
            window_average = np.mean(window)
              
            # Store the average of current window in moving average list
            moving_averages.append(window_average)

            # Shift window to right by one position
            i += 1
        
        for ind in range(window_size - 1):
            moving_averages.insert(0,np.nan)
            
        return moving_averages

    # Plot
    plt.plot(x,y,'o')
    
    #print(len(x),len(y))
    # Plot smooth function
    # bspl = splrep(x,y,k=3,s=5)
    # y_smooth = splev(x,bspl)
    # plt.plot(t,y_smooth)
    
    # y_av = movingaverage(y,5)
    # plt.plot(x, y_av)
    
    plt.ylim(0, 1.5)
    plt.xlim(0, 100)
    
    plt.xlabel('Distance from origin')
    plt.ylabel('Cell speed')
    plt.title('t={t}min'.format(t=t))
    
    
    
    #### 2D plot
    # X = positions_x
    # Y = positions_y

    # Z = total_speeds
    
    # plt.figure()
    
    # plt.tricontourf(X,Y,Z,levels=20)
    # plt.plot(X,Y,'ok')

    
    
    #### 3D plot
    
    # X = positions_x
    # Y = positions_y

    # Z = total_speeds
    
    # fig = plt.figure()
    
    # ax = fig.add_subplot(111, projection='3d')

    # ax.plot_trisurf(X,Y,Z,cmap='viridis',antialiased = True)

    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z');


    # plt.show()
    
    
if __name__ == '__main__':
    
    # Figure resolution
    mpl.rcParams['figure.dpi'] = 300
    
    # Plot style
    plt.style.use('ggplot')
    
    ribose = '50'
    
    # Data folder
    data_folder = '..\\output\\ribose' + ribose + '\\threads1\\output\\'
    
    # Save folder
    save_folder = '..\\results\\ribose' + ribose + '\\'    
    
    # Speed distribution histogram
    cell_speed_heatmap('output00000180',data_folder,save_folder)
    
    # files = os.listdir(data_folder)
    
    # for i in range(len(files)):
    #     if not re.search('_ECM\.mat', files[i]):
    #         continue
    #     cell_speed_heatmap(files[i].split('_')[0], data_folder,save_folder)
    
    # fig = plt.gcf()
    # fig.patch.set_alpha(0.0)

        