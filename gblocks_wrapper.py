"""
assuming that all alignment files end with ".aln"
"""

import os,sys

FILE_ENDING = ".aln"
seqtype = "--codon"

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python gblocks_wrapper.py inDIR outDIR"
		sys.exit()
	
	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	
	for i in os.listdir(inDIR):
		if i[-len(FILE_ENDING):] != FILE_ENDING: continue		
		cmd = "Gblocks "+inDIR+i+" -t=c"
		print cmd
		os.system(cmd)

