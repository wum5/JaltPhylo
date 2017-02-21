#!/usr/bin/python


import numpy, sys, os
from Bio import SeqIO


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "usage: python codemelScript.py inDIR codemlFile treeFile"
		sys.exit()
		
	inDIR = sys.argv[1]+"/"
	codeml = sys.argv[2]
	tree = sys.argv[3]
	
	for i in os.listdir(inDIR):
		if i[-3:] != "phy": continue		
		clusterID = i.split('.phy')[0]

		cmd = "mkdir "+inDIR+clusterID
		os.system(cmd)
		
		output = "codeml.ctl"

		with open(output,"w") as outfile:
			outfile.write('\tseqfile = '+i+'\t* sequence data filename\n')
			outfile.write('\ttreefile = '+tree+'\t* tree structure file name\n')
			outfile.write('\toutfile = '+clusterID+'_out\t* main result file name\n')
		
		open(output, "a").writelines([l for l in open(codeml).readlines()]) 
			
		outfile.close()


		cmd = "mv "+inDIR+i+" "+inDIR+clusterID
		os.system(cmd)
		cmd = "mv "+output+" "+inDIR+clusterID
		os.system(cmd)
		cmd = "cp "+tree+" "+inDIR+clusterID
		os.system(cmd)

