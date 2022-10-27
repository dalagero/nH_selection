#!/bin/bash

config='fit_config_nH_resid_flash_yasu.json'
filename='test'
setup='shape'

numVals_s22t13=1 #True number of values is this number+1
s22t13_low=0.05
s22t13_high=0.11

numVals_dm2ee=1 #True number of values is this number+1
dm2ee_low=0.0002
dm2ee_high=0.0055

#FIXME Add a thing here about which setup to do for the best fit
python fit.py $config --dm2ee 0.0025227 --avg-near --dm2ee-behavior 'free' --pulls all --shape --save-fit $filename

chmod 777 $filename.npz

for (( s22t13_counter=0; s22t13_counter<=$numVals_s22t13; s22t13_counter++ ))
do
    runVal_s22t13=$(echo "scale=10;$s22t13_low+$s22t13_counter*($s22t13_high-$s22t13_low)/$numVals_s22t13" | bc)
    sbatch -q debug -t 30 jobScript.sh $config $filename $setup $runVal_s22t13 $numVals_dm2ee $dm2ee_low $dm2ee_high
done

