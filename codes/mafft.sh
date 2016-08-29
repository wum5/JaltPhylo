#!/bin/bash

#PBS -N Jal-mafft
#PBS -l nodes=1:ppn=8,walltime=12:00:00,vmem=64gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load java
module load python

PATH=$PATH:/N/dc2/projects/jaltomt/softwares/mafft-7.222-without-extensions
cd /N/dc2/projects/jaltomt/scripts
OD=/N/dc2/projects/jaltomt/homologs/2nd/sub_folder_1


python mafft_wrapper.py $OD $OD .fa 4 DNA
