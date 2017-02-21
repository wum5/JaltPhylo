import numpy, sys, os
from Bio import SeqIO

FILE_ENDING = ".fa"
TRE_ENDING = ".tre"

def get_name(name):
	if 'JA' in name:
		return name.replace("_R_","")[:6]
	else:
		return name.replace("_R_","")[:5]

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "usage: python cluster_gene_ID.py inDIR treDIR outDIR"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	treDIR = sys.argv[2]+"/"
	outDIR = sys.argv[3]+"/"
	cluster, gene = [], []
	
	for i in os.listdir(treDIR):
		if i[-len(TRE_ENDING):] != TRE_ENDING: continue		
		handle = open(treDIR+i, "rU")
		tree = str(handle.readline())
		start = tree.find("Solyc")
		cluster.append(i.split('.')[0])
		gene.append(tree[start:start+18])
		handle.close()
	
	for i in os.listdir(inDIR):
		if i[-len(FILE_ENDING):] != FILE_ENDING: continue
		currID = i.split('.')[0]
		for j in xrange(len(gene)):
			if currID == cluster[j]:
				clusterID = gene[j]

		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")	
		
		sortedList = [f for f in sorted(infile, key=lambda x : x.id)]
		clusterID = clusterID+'.fa'

		check = []
		num = 0
		
		outfile = open(outDIR+clusterID,"w")
		for record in sortedList:
			outfile.write('>'+get_name(record.id)+'\n')
			outfile.write(str(record.seq)+'\n')
		
		outfile.close()
		handle.close()

		print clusterID
