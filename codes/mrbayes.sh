#!/bin/bash

#PBS -N Jal-mrbayes-20
#PBS -l nodes=1:ppn=1,walltime=24:00:00
#PBS -m bea
#PBS -M wum5@umail.iu.edu


PATH=$PATH:/N/dc2/projects/jaltomt/software/mrbayes_3.2.5/src/
export LD_LIBRARY_PATH=/N/dc2/projects/jaltomt/software/mrbayes_3.2.5/src/beagle-lib/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/N/dc2/projects/jaltomt/software/mrbayes_3.2.5/src/beagle-lib/lib/pkgconfig:$PKG_CONFIG_PATH


cd /N/dc2/projects/jaltomt/de_novo/e5_80/full_ortholog_0.2/MrBayes/sub_folder_20

for d in */ ; do
	cd $d
	mb *.nex
	cd ..
done