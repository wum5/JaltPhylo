"""
Read the concatenated fasta file and the mcl output
write individual fasta files for each cluster
"""

import sys,os
from Bio import SeqIO

MIN_LEN = 200

if __name__ =="__main__":
	if len(sys.argv) != 5:
		print "usage: write_fasta_files_from_mcl.py fasta mcl_outfile minimal_ingroup_taxa outDIR"
		sys.exit()
	
	fasta = sys.argv[1]
	mclfile = sys.argv[2]
	min_ingroup_taxa = int(sys.argv[3])
	outDIR = sys.argv[4]+"/"
	
	print "Reading mcl output file"
	clusterDICT = {} #key is seqID, value is clusterID
	count = 0
	with open(mclfile,"rU") as infile:
		for line in infile:
			if len(line) < 3: continue #ignore empty lines
			spls = line.strip().split('\t')
			
			#count number of ingroup taxa
			taxa = []
			for seqID in spls:
				if "@" in seqID: #only look at ingroup taxa
					taxonID = seqID[:6]
					if taxonID not in taxa:
						taxa.append(taxonID)
			if len(taxa) >= min_ingroup_taxa:
				count += 1
				clusterID = str(count)
				for seqID in spls:
					clusterDICT[seqID] = clusterID
					
	print "Reading the fasta file"
	handle = open(fasta,"rU")
	for record in SeqIO.parse(handle,"fasta"):
		seqid,seq = str(record.id),str(record.seq)
		if len(seq) >= MIN_LEN:
			try:
				clusterID = clusterDICT[seqid]
				with open(outDIR+"cc"+clusterID+".fa","a") as outfile:
					outfile.write(">"+seqid+"\n"+seq+"\n")
			except:
				pass
	handle.close()