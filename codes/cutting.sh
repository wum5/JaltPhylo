#!/bin/bash

#PBS -N Jal-cutting
#PBS -l nodes=1:ppn=4,walltime=8:00:00,vmem=240gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load java
module load python
cd /N/dc2/projects/jaltomt/softwares/phyutility
CD=/N/dc2/projects/jaltomt/homologs/1st/sub_folder_5
OD=/N/dc2/projects/jaltomt/homologs/2nd/sub_folder_5
PATH=$PATH:/N/dc2/projects/jaltomt/softwares/mafft-7.222-without-extensions:/N/dc2/projects/jaltomt/softwares/FastTree


python cut_long_branches_iter.py $CD $OD
