import sys, os


"""
find out unprocessed input files and save them into a new directory
raw files (inDir2) excluding the ones being processed (inDir1) 
"""

def find_unprocessed(ccID1,ccID2):
	ccID3 = list(set(ccID2) - set(ccID1))
	return ccID3

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "python cp_aln_files.py inDIR1 inDIR2 outDIR"
		sys.exit(0)
	
	inDIR1 = sys.argv[1]+"/"
	inDIR2 = sys.argv[2]+"/"
	outDIR = sys.argv[3]+"/"
	ccID1, ccID2, ccID3 = [], [], []
	for i in os.listdir(inDIR1): 
		if not i.endswith(".fa"): continue
		ccID1.append(i.split(".")[0])
	
	for i in os.listdir(inDIR2): 
		if not i.endswith(".fa"): continue
		ccID2.append(i.split(".")[0])
	
	ccID3 = find_unprocessed(ccID1,ccID2)

	for item in ccID3:
		file = inDIR2 + item
		os.system("cp "+file+".fa "+outDIR+"")
