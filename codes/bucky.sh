#!/bin/bash

#PBS -N Jal-bucky
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=64gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


PATH=$PATH:/N/dc2/projects/jaltomt/software/bucky-1.4.4/src
cd /N/dc2/projects/jaltomt/de_novo/e5_80/full_ortholog_0.2/MrBayes


for d in main_*; do cd $d; mbsum -n 501 -o mymb mymb.run?.t; cd ..; done

bucky -n 1000000 -sg -i jalt_infilelist
