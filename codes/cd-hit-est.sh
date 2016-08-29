#!/bin/bash

#PBS -N Jal-cdhit
#PBS -l nodes=1:ppn=4,walltime=8:00:00,vmem=32gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

module load python/2.7.9
module load cd-hit

CD=/N/dc2/projects/jaltomt/scripts
OD=/N/dc2/projects/jaltomt/trinity_outdir/cdhit_out

cp $CD/cd-hit-est.py $OD/
cd $OD

python cd-hit-est.py
