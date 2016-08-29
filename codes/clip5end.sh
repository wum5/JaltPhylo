#!/bin/bash

#PBS -N Jal-clipper
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=10gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


## The module is in Karst server
module use -a /N/soft/rhel6/modules/karst/LIFE-SCIENCES
module load fastx/0.0.13
cd /N/dc2/projects/jaltomt/rawdata
OD=/N/dc2/projects/jaltomt/rawdata/trimmed


for file in *shear1*; do fastx_trimmer -f 16 -m 50 -Q33 -i $file -o $OD/$file; done
