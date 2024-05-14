#!/bin/bash
#SBATCH --qos bbdefault
#SBATCH --ntasks 10
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes 1
#SBATCH --time 3000:0
#SBATCH --account spillf-systems-mechanobiology-health-disease

set -e

module purge
module load bluebear
module load Miniconda3/4.9.2
module load bear-apps/2022a
# module load FFmpeg

source activate python_imaging

cd python_imaging/

python parameter_analysis_simulations_plots.py

echo -e "Finished!"