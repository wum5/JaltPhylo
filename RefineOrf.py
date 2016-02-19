import os,sys
from Bio import SeqIO

'''
removing 3' and 5' ends of transcripts of Jaltomata sequences using 
the CDS boundary of tomato sequences on multiple sequence alignments
'''

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python RefineOrf.py inDir outDir"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	
	for i in os.listdir(inDIR):
		if i[-4:] != ".aln": continue
		clusterID = i.split(".")[0]
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")
	
		frame = False
		
		#for record in infile:
		record_dict = SeqIO.to_dict(infile)
		s = record_dict['Solyc'].seq
		gap = 0
		for j in range(0, len(s), 1):
			if s[j] == '-':
				gap += 1
			else:
				if gap%3 != 0:
					frame = True
				gap = 0
		
		if not frame:
			os.system("mv "+inDIR+i+" "+outDIR)
		

		handle.close()
