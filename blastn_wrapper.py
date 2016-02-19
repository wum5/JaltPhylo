"""
run blast search on all .fa files in the designated folder
change the name of database when needed
"""

import os,sys

if __name__ == "__main__":
	if len(sys.argv) != 1:
		print "usage: python blastn_wrapper.py"
		sys.exit()
	
	DIR = "./"
	for i in os.listdir(DIR):
		if i.endswith(".cd-hit"):
			cmd = "blastn -db all.cds.fa -query "+DIR+i+" -evalue 1 -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore' -out "+DIR+i[:-3]+"blastn -num_threads=11 -max_target_seqs 100"
			print cmd
			os.system(cmd)