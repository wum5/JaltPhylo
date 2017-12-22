#!/bin/bash

#PBS -N SNP_CALL
#PBS -l nodes=1:ppn=1,walltime=12:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load samtools
module load python/3.5.0

cd /N/dc2/projects/jaltomt/Phylogenomics/maptoref
OS=/N/dc2/projects/jaltomt/Softwares/mvftools
REF=ITAG2.4_genomic.fasta
REG=SL2.50ch01


samtools mpileup -uD -f $REF -r $REG JA0701_uniq.sort.bam JA0456_uniq.sort.bam JA0694_uniq.sort.bam \
JA0450_uniq.sort.bam JA0798_uniq.sort.bam JA0711_uniq.sort.bam JA0723_uniq.sort.bam JA0608_uniq.sort.bam \
JA0702_uniq.sort.bam JA0726_uniq.sort.bam JA0432_uniq.sort.bam JA0010_uniq.sort.bam JA0719_uniq.sort.bam \
JA0816_uniq.sort.bam | bcftools view -cg - > $REG'.vcf'

python3 $OS/vcf2mvf.py --vcf $REG'.vcf' --out $REG'.mvf' --lowdepth 10 --lowqual 30

python3 $OS/mvf_join.py --mvf SL2.50ch00.mvf SL2.50ch01.mvf SL2.50ch02.mvf SL2.50ch03.mvf \
SL2.50ch04.mvf SL2.50ch05.mvf SL2.50ch06.mvf SL2.50ch07.mvf SL2.50ch08.mvf SL2.50ch09.mvf \
SL2.50ch10.mvf SL2.50ch11.mvf SL2.50ch12.mvf --out combined.mvf
