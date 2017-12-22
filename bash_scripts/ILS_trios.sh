#!/bin/bash
#PBS -N Jal-trios
#PBS -l nodes=1:ppn=1,walltime=12:00:00,vmem=1gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

cd /N/dc2/projects/jaltomt/Phylogenomics/introgression
OD=/N/dc2/projects/jaltomt/Softwares/mvftools

## test ILS (which how many gene tree topologies (BABA+ABBA) are different from the species tree at particular internode)
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0702 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0726 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0432 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0010 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0719 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0816 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0726 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0432 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0010 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0719 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0816 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0726 JA0432 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0726 JA0010 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0726 JA0719 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0726 JA0816 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0432 JA0010 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0432 JA0719 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0432 JA0816 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0010 JA0719 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0010 JA0816 JA0723 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0719 JA0816 JA0723 --outgroup Solyc --windowsize 6201996

python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0010 JA0719 JA0816 --outgroup Solyc --windowsize 6201996

python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0702 JA0816 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0726 JA0816 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0432 JA0816 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0726 JA0816 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0432 JA0816 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0726 JA0432 JA0816 --outgroup Solyc --windowsize 6201996

python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0702 JA0719 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0726 JA0719 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0432 JA0719 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0726 JA0719 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0432 JA0719 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0726 JA0432 JA0719 --outgroup Solyc --windowsize 6201996

python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0702 JA0010 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0726 JA0010 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0608 JA0432 JA0010 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0726 JA0010 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0432 JA0010 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0726 JA0432 JA0010 --outgroup Solyc --windowsize 6201996

python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0726 JA0608 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0702 JA0432 JA0608 --outgroup Solyc --windowsize 6201996
python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0726 JA0432 JA0608 --outgroup Solyc --windowsize 6201996

python $OD/mvf_chromoplot.py --mvf transcriptome --samples JA0726 JA0432 JA0702 --outgroup Solyc --windowsize 6201996
