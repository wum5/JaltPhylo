#!/bin/bash

#PBS -N trinity
#PBS -l nodes=1:ppn=8,walltime=32:00:00,vmem=280gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load gcc/4.9.2
module load java/jre/1.8.0_73
module load fastqc
module load samtools/1.2
module load bowtie
module load trinityrnaseq


CD=/N/dc2/projects/jaltomt/rawdata/trimmed
OD=/N/dc2/projects/jaltomt/trinity_outdir/trinity_out
TG1=JA0456RP

Trinity --seqType fq --max_memory 80G --CPU 8 --full_cleanup --output $OD/$TG1'_trinity' --left $CD/$TG1'_shear1_p1.fastq' --right $CD/$TG1'_shear1_p2.fastq'

