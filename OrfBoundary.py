import os, sys
from Bio import SeqIO

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python OrfBoundary.py inDir outDir"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"

	nucleotides = ['A', 'T', 'C', 'G']
	
		
	for i in os.listdir(inDIR):
		if i[-3:] != ".fa": continue
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")
		codon_cov = []
		num = 0
	
		for record in infile:
			if num == 0:
				for x in range(0, len(record.seq)-2, 3):
					codon_cov.append(0)
			for x in range(len(codon_cov)):
				if record.seq[x*3] in nucleotides and record.seq[x*3+1] in nucleotides and record.seq[x*3+2] in nucleotides:
					codon_cov[x] += 1
				elif record.id == 'Solyc':
					codon_cov[x] -= 10
				else:
					codon_cov[x] = codon_cov[x]
			
			num += 1 
			
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")
						
		with open(outDIR+i,"w") as outfile:
			for key in infile:
				outfile.write('>' + key.id + '\n')
				seq = ''
				for x in range(len(codon_cov)):
					if codon_cov[x]>num/2:
						seq += str(key.seq[x*3:x*3+3])
				outfile.write(seq+'\n')

		handle.close()
		outfile.close()
		
		