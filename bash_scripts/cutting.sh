#!/bin/bash
#PBS -t 1-5
#PBS -N Jal-cutting
#PBS -l nodes=1:ppn=4,walltime=8:00:00,vmem=240gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load java
module load python
cd /N/dc2/projects/jaltomt/softwares/phyutility
CD=/N/dc2/projects/jaltomt/homologs/1st/sub_folder_$PBS_ARRAYID
OD=/N/dc2/projects/jaltomt/homologs/2nd/sub_folder_$PBS_ARRAYID
PATH=$PATH:/N/dc2/projects/jaltomt/softwares/mafft-7.222-without-extensions:/N/dc2/projects/jaltomt/softwares/FastTree


python cut_long_branches_iter.py $CD $OD
