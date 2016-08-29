########################################################
## mask pooly aligned regions using 15-bp sliding window
## with criteria: 1) > 3 mismatched within ingroups
##				  2) > 5 mismatched among (ingroups + outgroup1)
##				  3) > 7 mismatched among (ingroups + 2 outgroups) 
########################################################


import os, sys
from Bio import SeqIO


nucleotides = ['A', 'T', 'C', 'G']
dif = lambda seq1, seq2: sum(seq1[i] != seq2[i] for i in xrange(len(seq1)) if seq1[i] in nucleotides and seq2[i] in nucleotides)

def window_check(windows, **kw):
	maxEdit = []
	for x in xrange(len(windows)):
		num, allEdit = 0, 0.0
		for y in xrange(len(windows)):
			allEdit += dif(windows[x], windows[y])
			num += 1
		maxEdit.append(allEdit/num)
	return max(maxEdit)

def remove_by_mask(seq, mask_cutoff=0.2):
	mask = False
	masked_region = seq.count('N')
	masked_ratio = round(float(masked_region)/len(seq),4)
	if masked_ratio > mask_cutoff:
		mask = True
		print i+'\t'+str(masked_ratio*100)+'%\t(deleted)'
		return mask
	else:
		print i+'\t'+str(masked_ratio*100)+'%'

		
	
if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python SlidingWindows.py inDir outDir"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	win1_size, win2_size, win3_size = 15, 15, 15
		
	for i in os.listdir(inDIR):
		if i[-3:] != ".fa": continue
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")
		#windows1 is for all seqs, windows2 consider (ingroups+close outgroup), windows2 consider just ingroups
		windows1, windows2, windows3 = [], [], []
		num1, num2, num3, snp = 0, 0, 0, 0
	
		for record in infile:
			if num1 == 0:
				for x in xrange(len(record.seq)-win1_size):
					windows1.append([str(record.seq[x:x+win1_size])])
			else:
				for x in xrange(len(record.seq)-win1_size):
					windows1[x].append(str(record.seq[x:x+win1_size]))
						
			if 'Cap' not in str(record.id):
				if num2 == 0:
					for x in xrange(len(record.seq)-win2_size):
						windows2.append([str(record.seq[x:x+win2_size])])
				else:
					for x in xrange(len(record.seq)-win2_size):
						windows2[x].append(str(record.seq[x:x+win2_size]))
				num2 += 1 

			if 'Solyc' not in str(record.id) and 'Cap' not in str(record.id):
				if num3 == 0:
					for x in xrange(len(record.seq)-win3_size):
						windows3.append([str(record.seq[x:x+win3_size])])
				else:
					for x in xrange(len(record.seq)-win3_size):
						windows3[x].append(str(record.seq[x:x+win3_size]))
				num3 += 1 
											
			num1 += 1 

		mask_win1, mask_win2, mask_win3 = [], [], []
		for x in xrange(len(windows1)-win1_size):
			score = window_check(windows1[x])
			if score > 7:	# equal to >7 mutations in 15bps windows
				mask_win1.append(x)
		for x in xrange(len(windows2)-win2_size):
			score = window_check(windows2[x])
			if score > 5:	# equal to >5 mutations in 15bps windows
				mask_win2.append(x)
		for x in xrange(len(windows3)-win3_size):
			score = window_check(windows3[x])
			if score > 3:	# equal to >3 mutations in 15bps windows
				mask_win3.append(x)
			
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")
					
		with open(outDIR+i,"w") as outfile:
			for record in infile:
				outfile.write('>' + record.id + '\n')
				seq = record.seq
				for x in xrange(len(mask_win1)):
					start = mask_win1[x]
					seq = seq[:start]+'NNNNNNNNNNNNNNN'+seq[start+win1_size:]
				for x in xrange(len(mask_win2)):
					start = mask_win2[x]
					seq = seq[:start]+'NNNNNNNNNNNNNNN'+seq[start+win2_size:]
				for x in xrange(len(mask_win3)):
					start = mask_win3[x]
					seq = seq[:start]+'NNNNNNNNNNNNNNN'+seq[start+win3_size:]
				outfile.write(str(seq)+'\n')
		
		if remove_by_mask(seq) == True:
			os.remove(outDIR+i)
		
		handle.close()
		outfile.close()


