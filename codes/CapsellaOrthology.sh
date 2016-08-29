#!/bin/bash

#PBS -N Jal-Capsella
#PBS -l nodes=1:ppn=1,walltime=2:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load python
module load biopython

SF=/N/dc2/projects/jaltomt/scripts
OF=/N/dc2/projects/jaltomt/orthologs


python $SF/CapsellaOrtholog.py $OF/4th $OF/Tomato_Capsella.txt $OF/Capsicum.annuum.L_Zunla-1_v2.0_CDS.fa $OF/5th