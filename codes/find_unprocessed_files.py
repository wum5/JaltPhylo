################################################################
## find out unprocessed input files and save them into a new directory
## raw files (inDir2) excluding the ones being processed (inDir1) 
################################################################


import sys, os


def find_unprocessed(inDIR1,inDIR2):
	ccID1, ccID2, ccID3 = [], [], []		
	ccID1 = [i.split(".fa")[0] for i in os.listdir(inDIR1) if i.endswith(".fa")]
	ccID2 = [i.split(".fa")[0] for i in os.listdir(inDIR2) if i.endswith(".fa")]	
	ccID3 = list(set(ccID2) - set(ccID1))
	return ccID3


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "python find_unprocessed_files.py processedDIR originalDIR outDIR"
		sys.exit(0)
	
	inDIR1 = sys.argv[1]+"/"
	inDIR2 = sys.argv[2]+"/"
	outDIR = sys.argv[3]+"/"	

	ccID3 = find_unprocessed(inDIR1,inDIR2)

	for item in ccID3:
		file = inDIR2 + item
		os.system("cp "+file+".fa "+outDIR+"")
