"""
Read alignment after cutting
write individual alignment files for each ortholog

Since no taxon repeats
Shorten seq id to taxon id
"""

import sys,os,newick3,phylo3
from Bio import SeqIO

ORTHO_TREE_FILE_ENDING = ".tre"

def get_name(label):
	if '@' in label:
		return label.replace("_R_","")[:6]
	else:
		return label.replace("_R_","")[:5]

def get_front_labels(node):
	leaves = node.leaves()
	return [i.label for i in leaves]

if __name__ =="__main__":
	if len(sys.argv) != 4:
		print "usage: python write_alignments_from_orthologs.py alnDIR orthoTreeDIR outDIR"
		sys.exit()
	
	alnDIR = sys.argv[1]+"/"
	treDIR = sys.argv[2]+"/"
	outDIR = sys.argv[3]+"/"
	"""
	print "Reading the original fasta file"
	handle = open(fasta,"rU")
	#hash table of taxonID -> seqID -> seq
	seqDICT = {} #key is taxonID, value is seqID
	for seq_record in SeqIO.parse(handle,"fasta"):
		seqID,seq = str(seq_record.id),str(seq_record.seq)
		taxonID = get_name(seqID)
		if taxonID not in seqDICT:
			seqDICT[taxonID] = {} #key is taxonID, value is seq
			print "Adding sequence from",taxonID
		seqDICT[taxonID][seqID] = seq
	handle.close()
	"""
	for i in os.listdir(treDIR):
		if i[-len(ORTHO_TREE_FILE_ENDING):] == ORTHO_TREE_FILE_ENDING:
			clusterID= i.split(".")[0]
			spls = clusterID.split("_")
			#when using pep
			alnname = spls[0]+".fa.aln"
			#when using cds
			print i,alnname
			
			#read in the alignment into an dictionary
			seqDICT = {} #key is seqID, value is seq
			with open(alnDIR+alnname,"rU") as handle:
				for record in SeqIO.parse(handle,"fasta"):
					seqDICT[str(record.id)] = str(record.seq)

			#read in tree tips and write output alignment
			with open(treDIR+i,"r")as infile:
				intree = newick3.parse(infile.readline())
			labels = get_front_labels(intree)
			with open(outDIR+i.replace(ORTHO_TREE_FILE_ENDING,".aln"),"w") as outfile:
				for lab in labels:
					outfile.write(">"+get_name(lab)+"\n"+seqDICT[lab]+"\n")
