#!/bin/bash
#PBS -t 1-14
#PBS -N mrbayes
#PBS -l nodes=1:ppn=1,walltime=48:00:00
#PBS -m bea
#PBS -M wum5@umail.iu.edu


PATH=$PATH:/N/dc2/projects/jaltomt/softwares/mrbayes_3.2.5/src
export LD_LIBRARY_PATH=/N/dc2/projects/jaltomt/softwares/mrbayes_3.2.5/src/beagle-lib/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/N/dc2/projects/jaltomt/softwares/mrbayes_3.2.5/src/beagle-lib/lib/pkgconfig:$PKG_CONFIG_PATH


cd /N/dc2/projects/jaltomt/phylogeny/bucky_tre/sub_folder_$PBS_ARRAYID

for d in */ ; do
	cd $d
	mb *.nex
	cd ..
done