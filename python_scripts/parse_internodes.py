#!/usr/bin/env python
"""This script is to firstly filter gene trees (the concatenated set of the RAxML bipartition trees) 
by comparing the average bootstrap values against a cutoff threshold. Then, the script is to 
investigate the unsolved internode in my investiated speices tree,
python parse_internode.py -i genetrees.tre -b 80
"""

import argparse, re, os
import numpy as np
from Bio import Phylo
from cStringIO import StringIO


def bootstraps_extract(tree):
	rawlist = re.findall('\)\d+\:',tree)
	outlist = [int(x[1:][:-1]) for x in rawlist]
	return outlist


def filterTree(inDIR, cutoff1, cutoff2):
	outfile1 = open("trees_bootstrap_"+str(cutoff1), "w")
	outfile2 = open("genes_bootstrap_"+str(cutoff1), "w")
	for i in os.listdir(inDIR):
		if "RAxML_bipartitions" not in i: continue
		with open(inDIR+'/'+i, "r") as infile:
			tree = (infile.readline()).rstrip()
			bootstraps = bootstraps_extract(tree)
			if np.mean(bootstraps) >= cutoff1 and min(bootstraps) >= cutoff2:
				outfile1.write(tree+'\n')
				outfile2.write('.'.join(i.split('.')[1:])+'\n')
	outfile1.close()
	outfile2.close()
	
					
				
def main():
	parser = argparse.ArgumentParser(prog="parse_internode.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--inDIR', \
	help=("the directory having the RAxML bipartition trees"), required=True)
	parser.add_argument('-ab', '--average_bootstrap_cutoff', \
	help=("a designated average bootstrap cutoff"), default=50)
	parser.add_argument('-mb', '--minimum_bootstrap_cutoff', \
	help=("a designated minimum bootstrap cutoff"), default=10)
	args = parser.parse_args()
	
	inDIR = args.inDIR
	cutoff1 = int(args.average_bootstrap_cutoff)	
	cutoff2 = int(args.minimum_bootstrap_cutoff)	
	filterTree(inDIR, cutoff1, cutoff2)
	
	

if __name__ == "__main__":
	main()


