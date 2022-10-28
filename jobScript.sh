#!/bin/bash
#SBATCH -A dayabay 
#SBATCH -L SCRATCH 
#SBATCH -q regular 
#SBATCH -C haswell 
#SBATCH -t 01:00:00


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

if [ "$3" = "shape" ] || [ "$3" = "dm2" ]; then #only one s22t13 value (at a time for shape)
    echo "s22t13: $4"
    echo "number of dm2ee values: $5"
    echo "dm2ee_low: $6"
    echo "dm2ee_high: $7"

    for (( dm2ee_counter=0; dm2ee_counter<=$5; dm2ee_counter++ ))
    do
        runVal_dm2ee=$(echo "scale=10;$6+$dm2ee_counter*($7-$6)/$5" | bc)
        echo $runVal_dm2ee
        python make_chi2grid.py --config $1 --save $2 --set-up $3 --sin22theta $4 --dm2ee $runVal_dm2ee
    done

elif [ "$3" = "sin2" ] || [ "$3" = "rate" ]; then #only one dm2 value
    echo "dm2ee: $4"
    echo "number of sin2 values: $5"
    echo "sin2_low: $6"
    echo "sin2_high: $7"

    for (( sin2_counter=0; sin2_counter<=$5; sin2_counter++ ))
    do
        runVal_sin2=$(echo "scale=10;$6+$sin2_counter*($7-$6)/$5" | bc)
        echo $runVal_sin2
        python make_chi2grid.py --config $1 --save $2 --set-up $3 --sin22theta $runVal_sin2 --dm2ee $4
    done

else
    echo "No Bueno!! $3 is an incorrect input!!"
    
fi
