#!/bin/bash

#PBS -N Jal-paml
#PBS -l nodes=1:ppn=1,walltime=24:00:00
#PBS -m bea
#PBS -M wum5@umail.iu.edu

module load paml
cd /N/dc2/projects/jaltomt/paml/phylip_seqs/unprocessed


for file in Solyc*; do cd $file; codeml; cd ..; done
