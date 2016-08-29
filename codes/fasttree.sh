#!/bin/bash

#PBS -N Jal-fasttree
#PBS -l nodes=1:ppn=1,walltime=12:00:00,vmem=24gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load java
module load python
cd /N/dc2/projects/jaltomt/scripts

PATH=$PATH:/N/dc2/projects/jaltomt/softwares/FastTree
OD=/N/dc2/projects/jaltomt/homologs/1st/sub_folder_6


python fasttree_wrapper.py $OD $OD DNA
