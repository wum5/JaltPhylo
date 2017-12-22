#!/bin/bash

#PBS -N BUCKY
#PBS -l nodes=1:ppn=1,walltime=48:00:00,vmem=240gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

set -e
set -u
set -o pipefail
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/mrbayes-3.2.6/src
PATH=$PATH:/N/dc2/projects/jaltomt/Softwares/bucky-1.4.4/src


## pull out gene datasets with mean bootstrap value > 50
cd /N/dc2/projects/jaltomt/Phylogenomics/phylogeny/gene_tre
python ../../scripts/parse_internodes.py -i bipartition_tree -ab 50
cat genes_bootstrap_50 | while read line; do cp "input_seqs/"${line} bootstrap50_seqfiles/ ; done
python ../../scripts/seqformat_converter.py bootstrap50_seqfiles ../bucky_tre .fa .nex

## prepare Bucky input files
outgroup="Solyc"
cd /N/dc2/projects/jaltomt/Phylogenomics/phylogeny/bucky_tre/

text="\nbegin mrbayes;\n\toutgroup ${outgroup};\n\tset autoclose=yes nowarnings=yes;\n\tlset nst=6 rates=invgamma; [specifies a GTR+I+G]\n\tmcmc ngen=1000000 samplefreq=1000 printfreq=1000 savebrlens=yes filename=mymb;\n\tquit;\nend;"
for file in *.nex; do echo -e $text >> $file; done
python ../../scripts/directory_subpackage.py ./ 1 .nex

## run MrBayes (can be splited into 11 subset of genes for parallel and each takes ~64h)
num=$(ls -d sub_folder_* | wc -l)
for d in {801..1000}; do cd sub_folder_$d; mb *.nex; cd ..; done

## run BUCKy 
for d in sub_folder_{1..1190}; do echo ${d}"/mymb" >> bucky_infilelist; done
for d in sub_folder_{1..1190}; do cd $d; mbsum -n 501 -o mymb mymb.run?.t; cd ..; done
bucky -n 1000000 -sg -i bucky_infilelist
