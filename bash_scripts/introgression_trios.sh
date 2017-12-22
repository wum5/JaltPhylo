#!/bin/bash
#PBS -N Jal-trios
#PBS -l nodes=1:ppn=1,walltime=12:00:00,vmem=1gb
#PBS -m bea
#PBS -M wum5@umail.iu.edu

cd /N/dc2/projects/jaltomt/introgression

## test introgression
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0010 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0010 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0010 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0723 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0723 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0723 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0432 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0432 JA0456 --outgroup Solyc --windowsize 6201996 

python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0432 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0702 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0702 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0702 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0816 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0816 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0816 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0719 JA0694 --outgroup Solyc --windowsize 6201996 

python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0719 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0719 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0608 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0608 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0608 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0726 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0726 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0711 JA0726 JA0701 --outgroup Solyc --windowsize 6201996 

python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0010 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0010 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0010 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0723 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0723 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0723 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0432 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0432 JA0456 --outgroup Solyc --windowsize 6201996 

python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0432 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0702 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0702 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0702 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0816 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0816 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0816 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0719 JA0694 --outgroup Solyc --windowsize 6201996 

python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0719 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0719 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0608 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0608 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0608 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0726 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0726 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0798 JA0726 JA0701 --outgroup Solyc --windowsize 6201996 

python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0010 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0010 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0010 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0723 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0723 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0723 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0432 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0432 JA0456 --outgroup Solyc --windowsize 6201996 

python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0432 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0702 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0702 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0702 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0816 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0816 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0816 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0719 JA0694 --outgroup Solyc --windowsize 6201996 

python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0719 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0719 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0608 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0608 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0608 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0726 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0726 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0726 JA0701 --outgroup Solyc --windowsize 6201996 

python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0711 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0711 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0711 JA0701 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0798 JA0694 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0798 JA0456 --outgroup Solyc --windowsize 6201996 
python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples JA0450 JA0798 JA0701 --outgroup Solyc --windowsize 6201996 
