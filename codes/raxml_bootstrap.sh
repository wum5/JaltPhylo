#!/bin/bash

#PBS -N Jal-raxml_bootstrap
#PBS -l nodes=1:ppn=2,walltime=24:00:00,vmem=8gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

# if work on concatenated dataset
#cd /N/dc2/projects/jaltomt/phylogeny/concatenate
# if work on gene trees
cd /N/dc2/projects/jaltomt/phylogeny/gene_tre/sub_folder_4
PATH=$PATH:/N/dc2/projects/jaltomt/softwares/standard-RAxML-master


for file in *.phy; do raxmlHPC-PTHREADS-SSE3 -T 2 -f a -x 12345 -# 100 -p 12345 -m GTRCAT -o Solyc -s $file -n $file ; done
