import os,sys
from Bio import SeqIO
import numpy as np

'''
calculate the quantiles of alignment length and 
filter the data based on the quantile information
'''

def quantiles(size):
	q25 = np.percentile(size, 25)
	q50 = np.percentile(size, 50)
	q75 = np.percentile(size, 75)
	return q25, q50, q75	



if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: python SeqLength.py inDIR treDIR outDIR > outfile"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	#treDIR = sys.argv[2]+"/"
	#outDIR = sys.argv[3]+"/"
	
	size_array = []
	
	for i in os.listdir(inDIR):
		if i[-6:] != "aln-gb": continue
		handle = open(inDIR+i, "rU")
		records = list(SeqIO.parse(handle, "fasta"))
		handle.close()
		s = records[0].seq
		aln_size = len(s)
		size_array.append(aln_size)
		
	size = np.array(size_array)
	q25, q50, q75 = quantiles(size)
	

	print q75 

'''	
	for i in os.listdir(inDIR):
		if i[-6:] != "aln-gb": continue
		clusterID = i.split(".")[0]
		handle = open(inDIR+i, "rU")
		records = list(SeqIO.parse(handle, "fasta"))
		handle.close()
		s = records[0].seq
		aln_size = len(s)
		
		if aln_size > q75:
			print clusterID+'.nex'
			#os.system("cp "+inDIR+i+" "+outDIR)
			
			
		if aln_size > q75 and os.path.exists(treDIR+"RAxML_bestTree."+clusterID+".fa.aln-gb"):
			trefile = open(treDIR+"RAxML_bestTree."+clusterID+".fa.aln-gb", "rU")
			for line in trefile:
				line = line.rstrip()
				print line
'''
	