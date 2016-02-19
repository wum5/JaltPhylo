import sys
from Bio import SeqIO

"""
change the outgroup ID in the original fasta file
"""

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: python change_outgroupID.py inputfile > outputfile"
		sys.exit()
		
	inputfile = sys.argv[1]

	for seq_record in SeqIO.parse(inputfile, "fasta"):
		print '>'+seq_record.id.split(' ')[0]
		print seq_record.seq
