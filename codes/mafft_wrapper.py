"""
Takes a directory of fasta files

If there are >= 1000 sequences in the direction, use --auto
For fasta files with less than 1000 sequences, use the slower but much 
more accurate algorithm

Uncomment the com += "--anysymbol " line if there are "U" or any other unusual
charactors in the sequences
"""

import os,sys

if __name__ == "__main__":
	if len(sys.argv) != 6:
		print "usage: python mafft_wrapper.py inDIR outDIR infile_ending thread DNA/aa"
		sys.exit()
	
	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	file_end = sys.argv[3]
	thread = sys.argv[4]
	
	if sys.argv[5] == "aa":
		seqtype = "--amino"
	elif sys.argv[5] == "DNA":
		seqtype = "--nuc" 	
	else:
		print "Input data type: DNA or aa"
		sys.exit()
		
	done = os.listdir(outDIR)
	filecount = 0
	l = len(file_end)
	for i in os.listdir(inDIR):
		if i[-l:] == file_end and i+".aln" not in done:
			filecount += 1
			seqcount = 0 #record how many sequences in the fasta file
			maxlen = 0 #record the longest seq length
			seqlen = 0 #record the current seq length
			with open(inDIR+i,"r") as infile:
				for line in infile:
					if line[0] == ">":
						seqcount += 1
						maxlen = max(seqlen,maxlen)
						seqlen = 0
					else: seqlen += len((line.strip()).replace("-",""))
			if seqcount >= 1000 or maxlen >= 10000:
				alg = "--auto" #so that the run actually finishes!
			else: alg = "--genafpair --maxiterate 1000"
			com = "mafft "+alg+" "+seqtype+" --thread "+thread+" "
			#com += "--anysymbol " #use this when there are "U"s in aa sequences
			com += inDIR+i+" > "+outDIR+i+".aln"
			print com
			os.system(com)
	
	if filecount == 0:
		print "No file end with",file_end,"found"
				