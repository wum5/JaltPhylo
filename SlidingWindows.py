import os, sys
from Bio import SeqIO


nucleotides = ['A', 'T', 'C', 'G']


dif = lambda seq1, seq2: sum(seq1[i] != seq2[i] for i in range(len(seq1)) if seq1[i] in nucleotides and seq2[i] in nucleotides)

def window_check(windows, **kw):
	maxEdit = []
	for x in range(len(windows)):
		num, allEdit = 0, 0.0
		for y in range(len(windows)):
			allEdit += dif(windows[x], windows[y])
			num += 1
		maxEdit.append(allEdit/num)
	return max(maxEdit)
	
	
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python SlidingWindows.py inDir outDir"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"

	win1_size, win2_size = 15, 15
		
	for i in os.listdir(inDIR):
		if i[-3:] != ".fa": continue
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")
		windowsAll, windowsLess = [], []
		num, numI, snp = 0, 0, 0
	
		for record in infile:
			if num == 0:
				for x in range(len(record.seq)-win1_size):
					windowsAll.append([str(record.seq[x:x+win1_size])])
			else:
				for x in range(len(record.seq)-win1_size):
					windowsAll[x].append(str(record.seq[x:x+win1_size]))
						
			if 'Solyc' not in str(record.id):
				if numI == 0:
					for x in range(len(record.seq)-win2_size):
						windowsLess.append([str(record.seq[x:x+win2_size])])
				else:
					for x in range(len(record.seq)-win2_size):
						windowsLess[x].append(str(record.seq[x:x+win2_size]))
				numI += 1 
							
			num += 1 


		mask_win1, mask_win2 = [], []
		for x in range(len(windowsAll)-win1_size):
			score = window_check(windowsAll[x])
			if score > 3:	# equal to >5 mutations in 15bps windows
				mask_win1.append(x)
		
		for x in range(len(windowsLess)-win2_size):
			score = window_check(windowsLess[x])
			if score > 2:	# equal to >3 mutations in 15bps windows
				mask_win2.append(x)
			
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")
					
		with open(outDIR+i,"w") as outfile:
			for record in infile:
				outfile.write('>' + record.id + '\n')
				seq = record.seq
				for x in range(len(mask_win1)):
					start = mask_win1[x]
					seq = seq[:start]+'XXXXXXXXXXXXXXX'+seq[start+win1_size:]
				for x in range(len(mask_win2)):
					start = mask_win2[x]
					seq = seq[:start]+'XXXXXXXXXXXXXXX'+seq[start+win2_size:]
				outfile.write(str(seq)+'\n')

		masked_region = seq.count('X')
		masked_ratio = float(masked_region)/len(seq)
		print i+'\t'+str(masked_ratio*100)+'%'
		
		handle.close()
		outfile.close()
		
		if masked_ratio > 0.4:
			os.remove(outDIR+i)
			print i+' is removed'

