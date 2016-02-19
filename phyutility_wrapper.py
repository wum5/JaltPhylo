import sys,os
"""
remove columns with occupancy lower than MIN_COLUMN_OCCUPANCY
remove seqs shorter than MIN_CHR after filter columns
"""

MIN_CHR = 200

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python phyutility_wrapper.py inDIR outDIR MIN_COLUMN_OCCUPANCY DNA/aa"
		sys.exit(0)

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	MIN_COLUMN_OCCUPANCY = sys.argv[3]
	
	if sys.argv[4] == "aa":
		seqtype = " -aa"
	elif sys.argv[4] == "DNA":
		seqtype = "" 	
	else:
		print "Input data type: DNA or aa"
		sys.exit()

	for i in os.listdir(inDIR):
		if i[-3:] == "aln":
			cmd = "./phyutility"+seqtype+" -clean "+MIN_COLUMN_OCCUPANCY
			cmd += " -in "+inDIR+i+" -out "+outDIR+i+"-pht"
			print cmd
			os.system(cmd)
			
			#remove very short seqs
			seqid = ""
			infile = open(outDIR+i+"-pht","r")
			outfile = open(outDIR+i+"-cln","w")
			for line in infile:
				line = line.strip()
				if len(line) == 0: continue #skip empty lines
				if line[0] == ">":
					if seqid != "": #not at the first seqid
						if len(seq.replace("-","")) >= MIN_CHR:
							outfile.write(">"+seqid+"\n"+seq+"\n")
					seqid,seq = line.split(" ")[0][1:],""
				else: seq += line.strip()
			#process the last seq
			if len(seq.replace("-","")) >= MIN_CHR:
				outfile.write(">"+seqid+"\n"+seq+"\n")
			infile.close()
			outfile.close()
			
			#remove intermediate file
			os.system("rm "+outDIR+i+"-pht")
