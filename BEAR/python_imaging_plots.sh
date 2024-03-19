#!/bin/bash
#SBATCH --qos bbdefault
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks 3
#SBATCH --time 3000:0
#SBATCH --account spillf-systems-mechanobiology-health-disease

set -e

module purge
module load bluebear
module load Miniconda3/4.9.2
module load bear-apps/2022a

source activate python_imaging

cd user_projects/python_imaging/

# python simulations_plots.py
python parameter_analysis_simulations_plots.py

# python multi_cell_spheroid.py

# python multi_max_dist.py
# python multi_violin_plot.py

# python multi_spheroid_area.py
# python previous_current_position.py

echo -e "Finished!"