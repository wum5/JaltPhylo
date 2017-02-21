#!/bin/bash

#PBS -N Jal-bucky
#PBS -l nodes=1:ppn=1,walltime=40:00:00,vmem=80gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


PATH=$PATH:/N/dc2/projects/jaltomt/softwares/bucky-1.4.4/src
cd /N/dc2/projects/jaltomt/phylogeny/bucky_tre


#for d in Solyc*; do echo $d"/mymb" >> jalt_infilelist; done
#for d in Solyc*; do cd $d; mbsum -n 501 -o mymb mymb.run?.t; cd ..; done

bucky -n 1000000 -sg -i bucky_infilelist
