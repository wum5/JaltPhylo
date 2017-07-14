'''
Calculate the average bootstrap value for each gene tree
Prepare BUCKY inputs (with bootstraps; trees are selected based on average bootstrap cutoff)
'''

import sys, os, re


def avg_bootstrap(tre):
	mylist = re.findall(r'\[([^]]*)\]', tre)
	bt = [ int(x) for x in mylist ]
	mean = float(sum(bt))/len(bt)
	return mean


if __name__ == "__main__":
	if len(sys.argv) != 3:
			print "usage: python bucky_prepare.py inDir cutoff"
			sys.exit()

	inDIR = sys.argv[1]+"/"
	cutoff = int(sys.argv[2])
	tre_file = open("tre_file", "a")
	bs_file = open("bs_file", "a")
	bucky_infilelist = open("bucky_infilelist", "a")


	for i in os.listdir(inDIR):
		if 'bipartitionsBranchLabels' not in i: continue
		infile1 = open(inDIR+i, "r")
		ID = i.split('RAxML_bipartitionsBranchLabels.')[1]
		tre1 = infile1.readline()
		boot = avg_bootstrap(tre1)
		node_num = tre1.count('JA')+tre1.count('Solyc')
		if boot < cutoff or node_num != 15: continue

		infile2 = open(inDIR+'RAxML_bestTree.'+ID, "r")
		tre2 = infile2.readline()
		path1 = inDIR+'RAxML_bootstrap.'+ID+'\n'
		path2 = ID.split('.phy')[0]+'/mymb\n'

		tre_file.write(tre2)
		bs_file.write(path1)
		bucky_infilelist.write(path2)

		infile1.close()
		infile2.close()
	

	tre_file.close()
	bs_file.close()
	bucky_infilelist.close()
