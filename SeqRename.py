import numpy, sys, os
from Bio import SeqIO

FILE_ENDING = ".fa"

def get_name(name):
	if 'JA' in name:
		return name.replace("_R_","")[:6]
	else:
		return name.replace("_R_","")[:5]

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "usage: python rename.py inDIR outDIR IDfile"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	IDfile = open(sys.argv[3], "r")
	cluster, gene = [], []
	
	for line in IDfile:
		line = line.rstrip()
		cluster.append(line.split('\t')[0])
		try: gene.append(line.split('\t')[1])
		except: gene.append("None")
	
	for i in os.listdir(inDIR):
		if i[-len(FILE_ENDING):] != FILE_ENDING: continue
		currID = i.split('.')[0]
		for j in range(len(gene)):
			if currID == cluster[j]:
				clusterID = gene[j]

		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")	
		
		sortedList = [f for f in sorted(infile, key=lambda x : x.id)]
		clusterID = clusterID+'.fa'

		check = []
		num = 0
		
		outfile = open(outDIR+clusterID,"w")
		remove = False
		for record in sortedList:
			if num == 0:
				for j in range(len(record.seq)):
					check.append(1)
			for j in range(len(record.seq)):
				if record.seq[j] not in ["A", "T", "C", "G"]:
					check[j] = 0
				else: 
					''
				if len(str(record.seq)) < 200:
					remove = True
			num += 1

			outfile.write('>'+get_name(record.id)+'\n')
			outfile.write(str(record.seq)+'\n')
		
		outfile.close()
		handle.close()
	
		bases = int(check.count(1))
		if bases < 150 or remove:
			os.remove(outDIR+clusterID)
		
		print clusterID

		
