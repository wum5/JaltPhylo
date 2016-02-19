import sys, os
from Bio import SeqIO

"""
pull out all sequences for specified ID
"""

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python pick_ID_seq.py inDIR seqID > outfile"
		sys.exit()
		
	inDIR = sys.argv[1]+"/"
	seqID = sys.argv[2]

	
	for i in os.listdir(inDIR):
		for seq_record in SeqIO.parse(inDIR+i, "fasta"):
    			if seq_record.id == seqID:
				print '>'+i.split('.')[0]
    				print seq_record.seq