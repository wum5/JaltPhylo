#!/bin/bash

#PBS -N SNP_CALL
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=4gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load samtools
module load python/3.5.0

cd /N/dc2/projects/jaltomt/Phylogenomics/maptoref
OS=/N/dc2/projects/jaltomt/Softwares/mvftools
REF=ITAG2.4_genomic.fasta
REG=SL2.50ch01


#samtools mpileup -uD -f $REF -r $REG JA0701.sort.bam JA0456.sort.bam JA0694.sort.bam JA0450.sort.bam JA0798.sort.bam JA0711.sort.bam JA0723.sort.bam JA0608.sort.bam JA0702.sort.bam JA0726.sort.bam JA0432.sort.bam JA0010.sort.bam JA0719.sort.bam JA0816.sort.bam | bcftools view -cg - > $REG'.vcf'
python3 $OS/vcf2mvf.py --vcf $REG'.vcf' --out $REG'.mvf' --lowdepth 10 --lowqual 30
