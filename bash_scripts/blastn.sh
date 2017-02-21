#!/bin/bash

#PBS -N blastn
#PBS -l nodes=1:ppn=16,walltime=8:00:00,vmem=64GB
#PBS -m bea
#PBS -M wum5@umail.iu.edu

PATH=$PATH:/N/dc2/projects/jaltomt/softwares/blast/
OD=/N/dc2/projects/jaltomt/trinity_outdir/blast_out
SF=/N/dc2/projects/jaltomt/scripts

cp $SF/blastn_wrapper.py $OD/
cd $OD

cat *fa.cd-hit > all.cds.fa
makeblastdb -in all.cds.fa -parse_seqids -dbtype nucl -out all.cds.fa
python blastn_wrapper.py
