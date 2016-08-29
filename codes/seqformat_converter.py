from Bio import AlignIO
from Bio import Alphabet
import sys, os


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "usage: python seqformat_converter.py inDIR outDIR inputFormat(.fa/.phy)"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	FileEnding = sys.argv[3]

	for i in os.listdir(inDIR):
		if not i.endswith(FileEnding): continue
		clusterID = i.split(FileEnding)[0]
		if FileEnding == '.fa':
			AlignIO.convert(inDIR+i, "fasta", outDIR+clusterID+".phy", "phylip-sequential", alphabet=Alphabet.generic_dna)
		elif FileEnding == '.phy':
			AlignIO.convert(inDIR+i, "phylip-sequential", outDIR+clusterID+".fa", "fasta", alphabet=Alphabet.generic_dna)
		#AlignIO.convert(inDIR+i, "fasta", outDIR+clusterID+".nex", "nexus", alphabet=Alphabet.generic_dna)