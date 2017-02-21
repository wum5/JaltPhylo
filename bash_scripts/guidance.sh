#!/bin/bash
#PBS -t 1-50
#PBS -N guidance
#PBS -l nodes=1:ppn=4,walltime=30:00:00,vmem=4gb
#PBS -m abe


module load ruby
module load gcc/4.9.2
#module load perl
#module load bioperl

cd /N/dc2/projects/jaltomt/alignments/unprocessed/sub_folder_$PBS_ARRAYID
OD=/N/dc2/projects/jaltomt/alignments/temp
PATH=$PATH:/N/dc2/projects/jaltomt/softwares/mafft-7.222-without-extensions:/N/dc2/projects/jaltomt/softwares/prank-msa/src


for file in *fa; do mkdir $OD/$file; perl /N/dc2/projects/jaltomt/softwares/guidance.v2.0/www/Guidance/guidance.pl --seqFile $file --msaProgram PRANK --seqType codon --outDir $OD/$file --bootstraps 10 --proc_num 4; done
