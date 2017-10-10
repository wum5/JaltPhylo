#!/bin/bash

#PBS -N BUCKY
#PBS -l nodes=1:ppn=1,walltime=48:00:00,vmem=256gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

set -e
set -u
set -o pipefail
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/mrbayes-3.2.6/src
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/bucky-1.4.4/src

outgroup="Solyc"
cd /N/dc2/projects/jaltomt/phylogeny/bucky_tre/

## create input files for MrBayes
text="\nbegin mrbayes;\n\toutgroup ${outgroup};\n\tset autoclose=yes nowarnings=yes;\n\tlset nst=6 rates=invgamma; [specifies a GTR+I+G]\n\tmcmc ngen=1000000 samplefreq=1000 printfreq=1000 savebrlens=yes filename=mymb;\n\tquit;\nend;"
for file in *.nex; do echo -e $text >> $file; done
python directory_subpackage.py ./ 1 .nex

## run MrBayes (can be splited into subset of genes for parallel)
num=$(ls -d sub_folder_* | wc -l)
for d in sub_folder_{1..$num}; do cd $d; mb *.nex; cd ..; done

## run BUCKy 
for d in sub_folder{1..$num}; do echo $d"/mymb" >> bucky_infilelist; done
for d in sub_folder*; do cd $d; mbsum -n 501 -o mymb mymb.run?.t; cd ..; done
bucky -n 1000000 -sg -i bucky_infilelist

