# A hybrid computational model of cancer spheroid growth with ribose-induced collagen stiffening

This code was built upon the ECM framework developed by Metzcar et al. [[1](#references)] in PhysiCell [[2](#references)].


## Overview of custom folders and files
### File `simulations_loop.py`
Script to run multiple simulations automatically by giving the parameter values. The project (`proj`, e.g. `proj=ecm_density`) loaded must exist in the `user_projects` folder and the `data/proj` must also exist to save the simulations. The code creates a new folder for each simulation inside the `data/proj/` folder. Single simulations can also be initiated with this script. 

Alternatively, it is possible to manually load a project by using the command `make load PROJ={proj}` in the terminal and manually change the variables in the `config/PhysiCell_settings.xml` file. Then the project can be compiled using `make` and be executed using `./project`. By default, the simulation data will be saved in the `output` folder, but it can be changed in the `PhysiCell_settings.xml` file.


### Folder `user_projects`
This folder contains custom projects. In the current folder is present the project `ecm_density` which contains the necessary files to run the simulations:
- the `config` folder contains the `PhysiCell_settings.xml` file. This file contains the parameters used for the simulations.
- the `custom_modules` folder contains the `custom.h` and `custom.cpp` files. These files contain the custom functions used for the simulations of the model. This folder also contains `extracellular_matrix.h` and `extracellular_matrix.cpp` which contain the ECM element and ECM mesh class definitions and other initilization routine.

The files in this folder must be loaded to run the simulations. It is possible to do so by writing `make load PROJ={proj}` with `proj` being the name of the project, in this case `ecm_density`.


### Folder `data`
The data output of the simulations is saved in this folder when using the `simulations_loop.py` script. The current folder contains a subfolder named `ecm_density`. 


### Folder `results`
In this folder are saved the figures produced in the python scripts. The current folder contains a subfolder named `ecm_density`, which contains the subfolders `animations`, `images`, `plots` and `statistics`. 


### Folder `phython_imaging`
This folder contains the scripts to process the PhysiCell data and create figures. 

- `simulations_plots.py` contains the main code. It processes the simulation data in parallel and creates a DataFrame to analyse and plot the data.
- `simulation_data.py` creates and saves the DataFrame into a pickle file in the corresponding simulation folder in the `data/proj/` folder.
- `delaunay_function.py` and `spheroid_area_function.py` are called during the DataFrame creation. They compute the spheroid area and the Delaunay mean distance at each time point. They also produce figures that are saved in `results/proj/statistics`
- `box_plots.py`,  `heatmaps.py` and `plots_over_time.py` produce plot figures saved in the `results/proj/plots` folder.
- `cell_plus_environment_movie_maker.py` contains the function used to produce the simulation figures that show the cancer cell agents on top of the ECM density. The simulation figures are saved in the `results/proj/images` folder and can be used to produce a video, saved in the `results/proj/animations` folder.


## References
[1] Metzcar, J., Duggan, B.S., Fischer, B., Murphy, M., Heiland, R., Macklin, P, 2024. A simple framework for agent-based modeling with extracellular matrix. bioRxiv 2022.11.21.514608; https://doi.org/10.1101/2022.11.21.514608

[2] Ghaffarizadeh, A., Heiland, R., Friedman, S.H., Mumenthaler, S.M., and Macklin, P. PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellular Systems, PLoS Comput. Biol. 14(2): e1005991, 2018. https://dx.doi.org/10.1371/journal.pcbi.1005991





