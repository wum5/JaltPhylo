#!/bin/bash

#PBS -N Jal-astral
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

OS=/N/dc2/projects/jaltomt/softwares/ASTRAL-master
cd /N/dc2/projects/jaltomt/phylogeny/gene_tre
module load java

java -Xincgc -Xms3000m -Xmx3000m -XX:MaxPermSize=3000m -jar $OS/astral.4.10.9.jar -i tre_file -b bs_file -o astral_out.tre -r 100
