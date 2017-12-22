#!/bin/bash

#PBS -N Jal-astral
#PBS -l nodes=1:ppn=1,walltime=4:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

OS=/N/dc2/projects/jaltomt/Softwares/ASTRAL
module load java

cd /N/dc2/projects/jaltomt/Phylogenomics/phylogeny/gene_tre
# The folder contains all the best RAxML gene tree files for each gene
cat genes_bootstrap_50 | while read line; do cat "best_tree/RAxML_bestTree."${line} ; done \
>> ../astral_tre/genetrees.tre

# The boostrap folders contains all the bootstraps files generated for each gene trees
cat genes_bootstrap_50 | while read line; do echo \
"/N/dc2/projects/jaltomt/Phylogenomics/phylogeny/gene_tre/bootstrap/RAxML_bootstrap."${line} ; \
done >> ../astral_tre/bsfile.txt

cd /N/dc2/projects/jaltomt/Phylogenomics/phylogeny/astral_tre
java -Xincgc -Xms3000m -Xmx3000m -XX:MaxPermSize=3000m -jar $OS/astral.5.5.9.jar  \
-i genetrees.tre -b bsfile.txt -o astral_out.tre -r 100