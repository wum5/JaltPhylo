from Bio import AlignIO
from Bio import Alphabet
import sys, os

FileEnding = ".phy"

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python seqformat_converter.py inDIR outDIR"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	for i in os.listdir(inDIR):
		if not i.endswith(FileEnding): continue
		clusterID = i.split(FileEnding)[0]
		#AlignIO.convert(inDIR+i, "fasta", outDIR+clusterID+".phy", "phylip-sequential", alphabet=Alphabet.generic_dna)
		AlignIO.convert(inDIR+i, "phylip-sequential", outDIR+clusterID+".fa", "fasta", alphabet=Alphabet.generic_dna)
		#AlignIO.convert(inDIR+i, "fasta", outDIR+clusterID+".nex", "nexus", alphabet=Alphabet.generic_dna)