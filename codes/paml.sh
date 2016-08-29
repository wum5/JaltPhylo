#!/bin/bash

#PBS -N Jal-paml
#PBS -l nodes=1:ppn=1,walltime=8:00:00
#PBS -m bea
#PBS -M wum5@umail.iu.edu

module load paml
#cd /N/dc2/projects/jaltomt/de_novo/e5_80/updated_ortholog/post-guid/4th_preSWAMP/sub_folder_1
cd /N/dc2/projects/jaltomt/de_novo/e5_80/updated_ortholog/pranked_align/4th_preSWAMP/sub_folder_7


for file in Solyc*; do cd $file; codeml; cd ..; done
