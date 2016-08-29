"""
prepare multiple bash-script files for guidance alignment
"""

import os,sys

if __name__ =="__main__":
	if len(sys.argv) != 5:
		print "usage: python guidance_process.py DIR fileNum thread hour"
		sys.exit()
	
	DIR = sys.argv[1]+"/"
	num = int(sys.argv[2])
	thread = sys.argv[3]
	hour = sys.argv[4]

	for i in xrange(1, num+1):
		myfile = 'sub_folder_'+str(i)
		f = open(DIR+myfile+".sh", "w")
		f.write("#!/bin/bash\n")
		f.write("#PBS -N Jal-guidance-"+str(i)+"\n")
		f.write("#PBS -l nodes=1:ppn="+thread+",walltime="+hour+":00:00,vmem="+thread+"gb\n")
		f.write("#PBS -m bea\n")
		f.write("#PBS -M wum5@umail.iu.edu\n\n")

		f.write("module load ruby\n")
		f.write("module load gcc/4.9.2\n")
		f.write("module load perl\n")
		f.write("module load bioperl\n\n")

		f.write("cd /N/dc2/projects/jaltomt/alignments/unprocessed/"+myfile+"\n")
		f.write("OD=/N/dc2/projects/jaltomt/alignments/temp\n")
		f.write("PATH=$PATH:/N/dc2/projects/jaltomt/softwares/mafft-7.222-without-extensions:/N/dc2/projects/jaltomt/softwares/prank-msa/src\n\n")

		f.write("for file in *fa; do mkdir $OD/$file; perl /N/dc2/projects/jaltomt/softwares/guidance.v2.0/www/Guidance/guidance.pl --seqFile $file --msaProgram PRANK --seqType codon --outDir $OD/$file --bootstraps 10 --proc_num "+thread+"; done\n")

