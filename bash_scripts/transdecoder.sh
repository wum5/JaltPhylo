#!/bin/bash

#PBS -N transdecoder
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

module load transdecoder

cd /N/dc2/projects/jaltomt/trinity_outdir/trinity_out

for file in *fasta; do TransDecoder.LongOrfs -t $file; done
