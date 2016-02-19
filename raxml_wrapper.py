"""
Change the raxml command (RAXML_CMD) depend on which flavor of raxml and where the
executable is in your machine

Input: a dir of cleaned alignments in fasta format and end with "-cln"
Output: tree names are clusterID.raxml
"""

ALIGNMENT_FILE_ENDING = ".aln-cln"
#RAXML_CMD = "raxmlHPC-PTHREADS-SSE3"
RAXML_CMD = "raxmlHPC-PTHREADS-SSE3-icc"
import os,sys

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "python raxml_wrapper.py number_cores DNA/aa"
		sys.exit(0)
	
	num_cores = sys.argv[1]
	if sys.argv[2] == "aa": model = "PROTCATWAG"
	elif sys.argv[2] == "DNA": model = "GTRCAT" 	
	else:
		print "Input data type: DNA or aa"
		sys.exit()
	DIR = "./"
	
	l = len(ALIGNMENT_FILE_ENDING)
	for j in os.listdir(DIR):
		if j[-l:] != ALIGNMENT_FILE_ENDING or "RAxML_" in j:
			continue
		clusterID = j.split(".")[0]
		print clusterID
		cmd = RAXML_CMD+" -T "+num_cores+" -p 12345 -s "+j+" -n "+j+" -m "+model+" -o Capana Solyc"
		print cmd
		os.system(cmd)