#!/bin/bash

#PBS -N Jal-shear
#PBS -l nodes=1:ppn=1,walltime=12:00:00,vmem=20gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu


module load python/2.7.9

cd /N/dc2/projects/jaltomt/scripts
CD=/N/dc2/projects/jaltomt/rawdata/TEMPDIR
OD=/N/dc2/projects/jaltomt/rawdata

## require to change TG and BARCODE
TG=JA0798RP
BARCODE=CTGAAGCT

python shear.py --fq1 $CD/$TG'_p1.fq' --fq2 $CD/$TG'_p2.fq' --out1 $OD/$TG'_shear1_p1.fastq' --out2 $OD/$TG'_shear1_p2.fastq' --trimfixed 0:0 --trimqual 20:20  --filterlength 50 --trimpattern3 AGATC --trimpolyat 12 --trimambig --filterunpaired --filterambig 8 --filterlowinfo 0.5 --barcodes1 $BARCODE --execscythe /N/dc2/projects/jaltomt/softwares/scythe-master/scythe --platform TruSeq --filterfile $OD/$TG'_filtered' --tempdir $OD/TEMPDIR --log $OD/LOG --cleanup none

