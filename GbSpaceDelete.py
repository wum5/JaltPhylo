import sys, os
from Bio import SeqIO


if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python GbSpaceDelete.py inDIR outDIR"
		sys.exit()
		
	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	
	for i in os.listdir(inDIR):
		if i[-3:] != "-gb": continue
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")	
		
		sortedList = [f for f in sorted(infile, key=lambda x : x.id)]

		clusterID = i
		reMove = False
		with open(outDIR+clusterID,"w") as outfile:
			for record in sortedList:
				s = str(record.seq)
				s.replace(" ", "")
				if len(s) < 200: 
					reMove = True
				outfile.write('>'+str(record.id)+'\n')
				outfile.write(s+'\n')
			
		handle.close()
		outfile.close()

		if reMove:
			os.remove(outDIR+clusterID)