#!/bin/bash

#PBS -N Jal-phyutility
#PBS -l nodes=1:ppn=1,walltime=12:00:00,vmem=64gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load java
module load python

cd /N/dc2/projects/jaltomt/softwares/phyutility
OD=/N/dc2/projects/jaltomt/homologs/2nd/sub_folder_3

python phyutility_wrapper.py $OD $OD 0.1 DNA
