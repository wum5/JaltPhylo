import os, sys
from Bio import SeqIO


def transform(seq_file):
	myseq = {}
	for record in seq_file:
		myseq[str(record.id)] = str(record.seq)
	return myseq
	
def concatenate(seqs, seq_dict):
	for key in seq_dict:
		myseq = seq_dict[key]
		if len(myseq) < 200:
			break
		elif key in seqs.keys():
			seqs[key] += myseq
	return seqs	
		
		
		
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python ConcatSeqs.py inDir outDIR"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"


	chr = 'Solyc'
	transcriptome = {'JA0456':'','JA0450':'','JA0010':'','JA0608':'','JA0701':'','JA0694':'','JA0432':'','JA0702':'','JA0719':'','JA0723':'','JA0726':'','JA0816':'','JA0711':'','JA0798':'','Solyc':''}
	
	for i in range(1,13):
		seqs = {'JA0456':'','JA0450':'','JA0010':'','JA0608':'','JA0701':'','JA0694':'','JA0432':'','JA0702':'','JA0719':'','JA0723':'','JA0726':'','JA0816':'','JA0711':'','JA0798':'','Solyc':''}	
		gene_num = 0
		if i < 10:
			chr_num = chr+'0'+str(i)
		else:
			chr_num = chr+str(i)
		for j in os.listdir(inDIR):
			if j[-3:] != ".fa" or j[:7] != chr_num: continue
			handle = open(inDIR+j, "rU")
			infile = SeqIO.parse(handle, "fasta")
			myseq = transform(infile)
			seqs = concatenate(seqs, myseq)
			transcriptome = concatenate(transcriptome, myseq)
			gene_num += 1
			
		outfile = open(outDIR+"ch"+str(i)+".fa", "w")
		for key in sorted(seqs):
			outfile.write('>' + key + '\n')
			outfile.write(seqs[key]+'\n')
		outfile.close()
		
		print chr_num, gene_num
		
	outfile  = open(outDIR+"transcriptome"+".fa", "w")
	for key in sorted(transcriptome):
		outfile.write('>' + key + '\n')
		outfile.write(transcriptome[key]+'\n')
	outfile.close()

