#!/bin/bash

#PBS -N Jal-SlidingWindow
#PBS -l nodes=1:ppn=1,walltime=16:00:00,vmem=8gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load python
module load biopython

OF=/N/dc2/projects/jaltomt/alignments
SF=/N/dc2/projects/jaltomt/scripts/
python $SF/SlidingWindows.py $OF/2nd/sub_folder_3 $OF/3rd >> $OF/masked_statistic.txt
