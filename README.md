# A hybrid computational model of cancer spheroid growth with ribose-induced collagen stiffening

This code built upon the ECM framework developed by [[1](#references)] in PhysiCell


## Overview of custom folders and files

### Folder user_projects
This folder contains custom projects. In the current folder is present the project `ecm_density` which contains the necessary files to run the simulations:
-`config` folder, contains the `PhysiCell_settings.xml` file. This file contains the parameters used for the simulations.
-`custom_modules`, contains `custom.h` and `custom.cpp` files. These files contain the custom functions used for the simulations of the model. This folder also contains `extracellular_matrix.h` and `extracellular_matrix.cpp` which contain the ECM element and ECM mesh class definitions and other initilization routine.

The files in this folder must be loaded in order to run the simulations. It is possible to do so by writing `make load PROJ={proj}` with `proj` being the name of the project, in this case `ecm_density`.

### Folder data
In this folder is saved the data output of the simulations. The current folder contains a subfolder named `ecm_density`. 


### Folder results
In this folder are saved the figures produced in the python scripts. The current folder contains a subfolder named `ecm_density`, which contains the subfolders `animations`, `images`, `plots` and `statistics`. 


### Folder phython_imaging
This folder contains the scripts to process the PhysiCell data and create figures. 

-`simulations_plots.py` contains the main code. It processes the simulation data in parallel and creates a DataFrame to analyse and plot the data.
-`simulation_data.py` creates and saves the DataFrame into a pickle file in the corresponding simulation folder in the `data\proj\` folder.
 







### File simulations_loop.py
Script to run multiple simulations automatically by giving the parameter values. The project (`proj`, e.g. `proj=ecm_density`) we load must exist in the user_projects folder and a folder with the same name must exist in the data folder in order to save the simulations. The code creates a new folder for each simulation inside the data/`proj`/ folder.





