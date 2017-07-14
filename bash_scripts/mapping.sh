#!/bin/bash

#PBS -N Mapping
#PBS -l nodes=1:ppn=1,walltime=8:00:00,vmem=8gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


#module load bwa
module load gcc/4.9.2 
module load star/2.5.2b
module load samtools
module load python
module load java

cd /N/dc2/projects/jaltomt/Phylogenomics/maptoref
OD=/N/dc2/projects/jaltomt/Phylogenomics/rawdata
OS=/N/dc2/projects/jaltomt/Softwares/mvftools-master
REF=ITAG2.4_genomic.fasta
ThreadN=4

#bwa index $REF
#samtools faidx $REF 
#STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ITAG2.4_genomic.fasta --runThreadN $ThreadN

for SPECIES in JA0701 JA0456 JA0694 JA0450 JA0798 JA0711 JA0723 JA0608 JA0702 JA0726 JA0432 JA0010 JA0719 JA0816 
do 
	#bwa mem $REF $OD/$SPECIES'NR_shear1_p1.fastq' $OD/$SPECIES'NR_shear1_p2.fastq' -t $THRED_NUM > $SPECIES'NR.sam'
	#bwa mem $REF $OD/$SPECIES'RP_shear1_p1.fastq' $OD/$SPECIES'RP_shear1_p2.fastq' -t $THRED_NUM > $SPECIES'RP.sam'
	
	STAR --genomeDir ./ --readFilesIn $OD/$SPECIES'NR_shear1_p1.fastq' $OD/$SPECIES'NR_shear1_p2.fastq' --outFileNamePrefix $SPECIES'NR.' --outSAMorder PairedKeepInputOrder --outReadsUnmapped Fastx --outSAMtype SAM --runThreadN $ThreadN
	STAR --genomeDir ./ --readFilesIn $OD/$SPECIES'RP_shear1_p1.fastq' $OD/$SPECIES'RP_shear1_p2.fastq' --outFileNamePrefix $SPECIES'RP.' --outSAMorder PairedKeepInputOrder --outReadsUnmapped Fastx --outSAMtype SAM --runThreadN $ThreadN

	samtools view -bS $SPECIES'NR.Aligned.out.sam' > $SPECIES'NR.bam'
	samtools view -bS $SPECIES'RP.Aligned.out.sam' > $SPECIES'RP.bam'  

	samtools sort -@ $ThreadN -m 5000000000 $SPECIES'NR.bam' $SPECIES'NR.sort'
	samtools sort -@ $ThreadN -m 5000000000 $SPECIES'RP.bam' $SPECIES'RP.sort'

	samtools merge $SPECIES'.sort.bam' $SPECIES'NR.sort.bam' $SPECIES'RP.sort.bam'
	#java -Xmx4G -jar /N/soft/rhel6/picard/picard-tools-1.107/MarkDuplicates.jar I=$SPECIES'.sort.bam' O=$SPECIES'.sort_mdup.bam' M=$SPECIES'_MarkDuplicate.txt' REMOVE_DUPLICATES=false AS=true
	samtools index $SPECIES'.sort.bam'
	
done

