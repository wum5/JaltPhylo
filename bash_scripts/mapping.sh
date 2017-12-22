#!/bin/bash

#PBS -N Mapping
#PBS -l nodes=1:ppn=4,walltime=12:00:00,vmem=16gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load star
module load samtools
module load python
module load java

cd /N/dc2/projects/jaltomt/Phylogenomics/maptoref
OD=/N/dc2/projects/jaltomt/Phylogenomics/KEY_FILE/rawdata
OS=/N/dc2/projects/jaltomt/Softwares/mvftools
REF=ITAG2.4_genomic.fasta
ThreadN=4

samtools faidx $REF 
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ITAG2.4_genomic.fasta --runThreadN $ThreadN

for SPECIES in JA0701 JA0456 JA0694 JA0450 JA0798 JA0711 JA0723 JA0608 JA0702 JA0726 JA0432 JA0010 JA0719 JA0816
do 
	
	STAR --genomeDir ./ --readFilesIn $OD/${SPECIES}'NR_shear1_p1.fastq' $OD/${SPECIES}'NR_shear1_p2.fastq' \
	--outFileNamePrefix ${SPECIES}'NR.' --outSAMorder PairedKeepInputOrder --outReadsUnmapped Fastx \
	--outSAMtype SAM --runThreadN $ThreadN
	STAR --genomeDir ./ --readFilesIn $OD/${SPECIES}'RP_shear1_p1.fastq' $OD/${SPECIES}'RP_shear1_p2.fastq' \
	--outFileNamePrefix ${SPECIES}'RP.' --outSAMorder PairedKeepInputOrder --outReadsUnmapped Fastx \
	--outSAMtype SAM --runThreadN $ThreadN

done

## here we only extract perfectly uniquely mapped reads (with flag "-q 255"; for STAR generated SAM file)
for SPECIES in JA0701 JA0456 JA0694 JA0450 JA0798 JA0711 JA0723 JA0608 JA0702 JA0726 JA0432 JA0010 JA0719 JA0816
do 

	samtools view -bS -q 255 ${SPECIES}'NR.Aligned.out.sam' | samtools sort -@ $ThreadN -m 5000000000 > ${SPECIES}'NR_uniq.sort.bam'
	samtools view -bS -q 255 ${SPECIES}'RP.Aligned.out.sam' | samtools sort -@ $ThreadN -m 5000000000 > ${SPECIES}'RP_uniq.sort.bam'

	samtools merge ${SPECIES}'_uniq.sort.bam' ${SPECIES}'NR_uniq.sort.bam' ${SPECIES}'RP_uniq.sort.bam'
	samtools index ${SPECIES}'_uniq.sort.bam'
	
done

