from Bio import AlignIO
from Bio import Alphabet
import sys, os


if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "usage: python seqformat_converter.py inDIR outDIR inputFormat(.fa/.phy/.nex) outputFormat(.fa/.phy/.nex)"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	inputFormat = sys.argv[3]
	outputFormat = sys.argv[4]

	for i in os.listdir(inDIR):
		if not i.endswith(inputFormat): continue
		clusterID = i.split(inputFormat)[0]
		if inputFormat == '.fa' and outputFormat == '.phy':
			AlignIO.convert(inDIR+i, "fasta", outDIR+clusterID+".phy", "phylip-sequential", alphabet=Alphabet.generic_dna)
		elif inputFormat == '.phy' and outputFormat == '.fa':
			AlignIO.convert(inDIR+i, "phylip-sequential", outDIR+clusterID+".fa", "fasta", alphabet=Alphabet.generic_dna)
		elif inputFormat == '.phy' and outputFormat == '.nex':
			AlignIO.convert(inDIR+i, "phylip-sequential", outDIR+clusterID+".nex", "nexus", alphabet=Alphabet.generic_dna)

