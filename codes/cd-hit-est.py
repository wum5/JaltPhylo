"""
run blast search on all .fa files in the designated folder
change the name of database when needed
"""

import os,sys

if __name__ == "__main__":
	if len(sys.argv) != 1:
		print "usage: python cd-hit-est.py"
		sys.exit()
	
	DIR = "./"
	for i in os.listdir(DIR):
		if i.endswith(".fasta") or i.endswith(".fa"):
			cmd = "cd-hit-est -i "+DIR+i+" -o "+DIR+i+".cd-hit -c 0.99 -n 10 -r 0 -T 4"
			print cmd
			os.system(cmd)
