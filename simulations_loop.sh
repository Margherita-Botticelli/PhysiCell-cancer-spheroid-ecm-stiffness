#!/bin/bash
#SBATCH --qos bbdefault
#SBATCH --ntasks 50
#SBATCH --time 3000:0
#SBATCH --account=spillf-systems-mechanobiology-health-disease

set -e

module purge; module load bluebear
module load GCC/9.3.0
module load Miniconda3/4.9.2
module load bear-apps/2022a

source activate python_imaging

python ./simulations_loop.py

echo "Finished!"