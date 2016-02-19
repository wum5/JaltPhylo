import numpy, sys, os
from Bio import SeqIO

def findStopCodons(orf, num):
	stopCodon = False
	for i in range(0, len(orf)-3, 3):
		codon = orf[i:i + 3]
		if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
			stopCodon = True
			num += 1
			break
	return stopCodon, num

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python StopCodon.py inDIR FileENding"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	FileENding = sys.argv[2]
	
	
	for i in os.listdir(inDIR):
		if i[-len(FileENding):] != FileENding: continue
		clusterID = i.split(".")[0]
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")	
		
		check = False
		num = 0
		
		for record in infile:
			orf = str(record.seq)
			check_new, num = findStopCodons(orf, num)
			check = max(check, check_new)
		
		if check:
			print clusterID, num
	
		handle.close()