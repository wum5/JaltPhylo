#!/bin/bash

#PBS -N Jal-FQC
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=8gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

module load java
module load fastqc

CD=/N/dc2/projects/jaltomt/rawdata/trimmed
OD=/N/dc2/projects/jaltomt/scripts
TG=JA0608RP

fastqc -o $OD -t 1 --nogroup $CD/$TG'_shear1_p1.fastq'

