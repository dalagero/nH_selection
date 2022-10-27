#!/bin/bash
#SBATCH -A dayabay 
#SBATCH -L SCRATCH 
#SBATCH -q regular 
#SBATCH -C haswell 
#SBATCH -t 00:30:00


cd /global/homes/d/dalagero/dyb-event-selection/

module load root
module load python
pip install --user -e .
export LBNL_FIT_P17B="1"
export OMP_NUM_THREADS=1
cd dyb_analysis/fitter


echo "config file: $1"
echo "filename: $2"
echo "setup: $3"
echo "s22t13: $4"
echo "dm2ee: $5"

python make_chi2grid.py --config $1 --save $2 --set-up $3 --sin22theta $4 --dm2ee $5
