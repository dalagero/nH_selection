#!/bin/bash -l
source /afs/ihep.ac.cn/users/b/beda/setup_nvwa.sh
cd /dybfs/users/odalager/nH
root -b -q "wrapper_ibds2000.c($1)"

