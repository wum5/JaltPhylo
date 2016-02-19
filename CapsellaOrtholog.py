import sys, os
from Bio import SeqIO

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "usage: python CapsellaOrtholog.py inDIR IDfile fastaFile outDIR"
		sys.exit()
		
	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[4]+"/"
	IDfile = open(sys.argv[2], "r")
	all, Capsella = [], []
	
	for line in IDfile:
		line = line.rstrip()
		all.append(line)
		Capsella.append(line.split('\t')[0])
	
	for i in os.listdir(inDIR):
		capID = ''
		if i[-3:] != "-gb": continue
		geneID = i.split('.masked')[0]
		for x in range(len(all)):
			if geneID in all[x]:
				capID = Capsella[x]
				print capID
				break
		
		if capID == '': continue
		
		print capID
		
		fastaFile = open(sys.argv[3], "r")
		handle = SeqIO.parse(fastaFile, "fasta")	
		for record in handle:
			currID = str(record.id)
			if capID in currID:
				capSeq = str(record.seq)
				break
		
		with file(inDIR+i, 'r') as original: data = original.read()
		with file(outDIR+i, 'w') as modified: modified.write('>Capana\n' + capSeq + '\n' + data)
