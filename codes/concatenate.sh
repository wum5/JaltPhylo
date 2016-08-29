#!/bin/bash

#PBS -N Jal-concatenate-masked
#PBS -l nodes=1:ppn=1,walltime=5:00:00,vmem=12gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load java
module load python
module load biopython
cd /N/dc2/projects/jaltomt/software/phyutility


python concatenate_matrices.py /N/dc2/projects/jaltomt/de_novo/e5_80/updated_ortholog/pranked_align/post2_SWAMP/GapRemoved 200 11 dna /N/dc2/projects/jaltomt/de_novo/e5_80/updated_ortholog/PAML/jalt_masked_6000NEW
