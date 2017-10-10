#!/bin/bash

#PBS -N Jal-astral
#PBS -l nodes=1:ppn=1,walltime=96:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

OS=/N/dc2/projects/jaltomt/Softwares/ASTRAL
cd /N/dc2/projects/jaltomt/Phylogenomics/phylogeny/gene_tre
module load java

# The folder contains all the best RAxML gene tree files for each gene
find ./ -name "*bestTree*" | sort -V | xargs cat > ../genetrees.tre
cd ..
# The boostrap folders contains all the bootstraps files generated for each gene trees
find bootstrap/ -name "*bootstrap*" | sort -V > bsfile.txt

java -Xincgc -Xms3000m -Xmx3000m -XX:MaxPermSize=3000m -jar $OS/astral.4.10.12.jar  \
-i genetrees.tre -b bsfile.txt -o astral_out.tre -r 100