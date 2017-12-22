#!/bin/bash

#PBS -N Jal-raxml_bootstrap
#PBS -l nodes=1:ppn=8,walltime=24:00:00,vmem=32gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load raxml/gnu/8.2.11 

cd /N/dc2/projects/jaltomt/Phylogenomics/phylogeny/gene_tre/input_seqs
for file in *.fa; do raxmlHPC-PTHREADS -T 8 -f a -x 12345 -# 100 -p 12345 \
-m GTRGAMMA -o Solyc -s ${file} -n ${file}; done

mv RAxML_bipartitions.* ../bipartition_tree
mv RAxML_bestTree.* ../best_tree
mv RAxML_bootstrap.* ../boostrap
mv RAxML* ../other_raxmlout