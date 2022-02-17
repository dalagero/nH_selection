#!/bin/bash -l
source /afs/ihep.ac.cn/users/b/beda/setup_nvwa.sh
cd /dybfs/users/odalager/nH
root -b -q "wrapper_singles400.c($1)"

