#!/bin/bash

#PBS -N MS_SIM
#PBS -l nodes=1:ppn=1,walltime=24:00:00,vmem=8gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

OD=/N/dc2/projects/jaltomt/Softwares/msdir 
cd /N/dc2/projects/jaltomt/Phylogenomics/simulation


# Each coalescent unit is 1 Mya since 4Ne = 4 × 10^5 individuals × 2.5 generations per year
$OD/ms 14 1000000000 -s 1 -I 14 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \
-ej 0.6852423 13 12 -ej 0.8062568 14 12 -ej 0.7957158 11 10 \
-ej 0.913646 10 9 -ej 0.9747881 9 8 -ej 1.052372 12 8 \
-ej 1.154937 8 7 -ej 0.7467259 6 5 -ej 1.287527 7 5 \
-ej 1.583158 5 4 -ej 1.524519 2 1 -ej 1.687863 3 1 \
-ej 2.901874 4 1 > ms_sim_Ne1e5.txt

# Each coalescent unit is 2 Mya by assuming Ne = 2 × 10^5 individuals 
$OD/ms 14 1000000000 -s 1 -I 14 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \
-ej 0.34262115 13 12 -ej 0.4031284 14 12 -ej 0.3978579 11 10 \
-ej 0.456823 10 9 -ej 0.48739405 9 8 -ej 0.526186 12 8 \
-ej 0.5774685 8 7 -ej 0.37336295 6 5 -ej 0.6437635 7 5 \
-ej 0.791579 5 4 -ej 0.7622595 2 1 -ej 0.8439315 3 1 \
-ej 1.450937 4 1 > ms_sim_Ne2e5.txt

# Each coalescent unit is 4 Mya by assuming Ne = 4 × 10^5 individuals 
$OD/ms 14 1000000000 -s 1 -I 14 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \
-ej 0.171310575 13 12 -ej 0.2015642 14 12 -ej 0.19892895 11 10 \
-ej 0.2284115 10 9 -ej 0.243697025 9 8 -ej 0.263093 12 8 \
-ej 0.28873425 8 7 -ej 0.186681475 6 5 -ej 0.32188175 7 5 \
-ej 0.3957895 5 4 -ej 0.38112975 2 1 -ej 0.42196575 3 1 \
-ej 0.7254685 4 1 > ms_sim_Ne4e5.txt

# Each coalescent unit is 0.5 Mya by assuming Ne = 0.5 × 10^5 individuals 
$OD/ms 14 1000000000 -s 1 -I 14 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \
-ej 1.3704846 13 12 -ej 1.6125136 14 12 -ej 1.5914316 11 10 \
-ej 1.827292 10 9 -ej 1.9495762 9 8 -ej 2.104744 12 8 \
-ej 2.309874 8 7 -ej 1.4934518 6 5 -ej 2.575054 7 5 \
-ej 3.166316 5 4 -ej 3.049038 2 1 -ej 3.375726 3 1 \
-ej 5.803748 4 1 > ms_sim_Ne5e4.txt

# Each coalescent unit is 0.25 Mya by assuming Ne = 0.25 × 10^5 individuals 
$OD/ms 14 1000000000 -s 1 -I 14 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \
-ej 2.7409692 13 12 -ej 3.2250272 14 12 -ej 3.1828632 11 10 \
-ej 3.654584 10 9 -ej 3.8991524 9 8 -ej 4.209488 12 8 \
-ej 4.619748 8 7 -ej 2.9869036 6 5 -ej 5.150108 7 5 \
-ej 6.332632 5 4 -ej 6.098076 2 1 -ej 6.751452 3 1 \
-ej 11.607496 4 1 > ms_sim_Ne25e3.txt


# I observed 31 non-synonymous mutations associated with the patter
# 0000110101111
python ../scripts/phyloGWAS_pval.py -i Jalt_noSolyc_codon \
-m ms_sim_Ne1e5.txt -p 0000110101111 -n 31 > output_Ne1e5.txt

python ../scripts/phyloGWAS_pval.py -i Jalt_noSolyc_codon \
-m ms_sim_Ne2e5.txt -p 0000110101111 -n 31 > output_Ne2e5.txt

python ../scripts/phyloGWAS_pval.py -i Jalt_noSolyc_codon \
-m ms_sim_Ne4e5.txt -p 0000110101111 -n 31 > output_Ne4e5.txt

python ../scripts/phyloGWAS_pval.py -i Jalt_noSolyc_codon \
-m ms_sim_Ne5e4.txt -p 0000110101111 -n 31 > output_Ne5e4.txt

python ../scripts/phyloGWAS_pval.py -i Jalt_noSolyc_codon \
-m ms_sim_Ne25e3.txt -p 0000110101111 -n 31 > output_Ne25e3.txt

