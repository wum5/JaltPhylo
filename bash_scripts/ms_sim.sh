#!/bin/bash

#PBS -N MS_SIM
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

OD=/N/dc2/projects/jaltomt/Softwares/msdir 
cd /N/dc2/projects/jaltomt/Phylogenomics/simulation

$OD/ms 14 1000000000 -s 1 -I 14 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \
-ej 0.687437 13 12 -ej 0.8088709 14 12 -ej 0.7992298 11 10 \
-ej 0.9178291 10 9 -ej 0.9788585 9 8 -ej 1.056161 12 8 \
-ej 1.158889 8 7 -ej 0.7485802 6 5 -ej 1.291403 7 5 \
-ej 1.587324 5 4 -ej 1.527497 2 1 -ej 1.690996 3 1 \
-ej 2.906452 4 1 > ms_sim.txt

# I observed 31 non-synonymous mutations associated with the patter
# 0000110101111
python ../scripts/phyloGWAS_pval.py -i Jalt_noSolyc_codon \
-m ms_sim.txt -p 0000110101111 -n 31 > output.txt
