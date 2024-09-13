#!/bin/bash
#SBATCH --qos bbdefault
#SBATCH --ntasks 20
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes 1
#SBATCH --time 3000:0
#SBATCH --account spillf-systems-mechanobiology-health-disease

set -e

module purge
module load bluebear
module load bear-apps/2022b
module load Miniforge3/24.1.2-0

source activate python_imaging

cd python_imaging/

python simulation_plots.py

echo -e "Finished!"