#!/bin/bash

#PBS -N Jal-raxml
#PBS -l nodes=1:ppn=4,walltime=24:00:00,vmem=24gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


OD=/N/dc2/projects/jaltomt/alignments/4th
PATH=$PATH:/N/dc2/projects/jaltomt/softwares/standard-RAxML-master

cd /N/dc2/projects/jaltomt/scripts
cp raxml_wrapper.py $OD
cd $OD


python raxml_wrapper.py 4 DNA
