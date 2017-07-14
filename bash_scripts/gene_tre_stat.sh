#!/bin/bash

#PBS -N gene_tree_stat
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=4gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

module load python
module load biopython

CD=/N/dc2/projects/jaltomt/Phylogenomics/scripts
cd /N/dc2/projects/jaltomt/Phylogenomics/phylogeny/gene_tre

python $CD/compare_topology.py genesTree.txt speciesTree.txt