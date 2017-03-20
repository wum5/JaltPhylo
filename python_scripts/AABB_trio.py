import itertools, sys


if __name__ == "__main__":
	if len(sys.argv) != 1:
		print "usage: python AABB_trios.py"
		sys.exit()


	green = ['JA0711','JA0798']
	orange = ['JA0010', 'JA0723','JA0432','JA0702','JA0816','JA0719', 'JA0608', 'JA0726']
	black = ['JA0694','JA0456','JA0701']
	red = ['JA0450']

	trios1 = list(itertools.product(green,orange,black))
	trios2 = list(itertools.product(red,orange,black))
	trios3 = list(itertools.product(red,green,black))
	trios = trios1+trios2+trios3


	f = open("trios.sh", "w")
	f.write("#!/bin/bash\n")
	f.write("#PBS -N Jal-trios\n")
	f.write("#PBS -l nodes=1:ppn=1,walltime=12:00:00,vmem=4gb\n")
	f.write("#PBS -m bea\n")
	f.write("#PBS -M wum5@umail.iu.edu\n\n")

	f.write("cd /N/dc2/projects/jaltomt/introgression\n\n")


	pre = 'python /N/dc2/projects/jaltomt/softwares/mvftools/mvf_chromoplot.py --mvf transcriptome --samples '
	pos = '--outgroup Solyc --windowsize 6204538'


	for key in trios:
		curr = ''
		for x in key:
			curr += x + ' '
		f.write(pre+curr+pos+"\n")

