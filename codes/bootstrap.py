import sys, os, re

def avg_bootstrap(tre):
	mylist = re.findall(r'\[([^]]*)\]', tre)
	bt = [ int(x) for x in mylist ]
	mean = float(sum(bt))/len(bt)
	return mean


'''
calculate the average bootstrap value for each gene tree
'''

if __name__ == "__main__":
	if len(sys.argv) != 4:
			print "usage: python bootstrap.py inDir outDIR cutoff"
			sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"
	cutoff = int(sys.argv[3])

	for i in os.listdir(inDIR):
		if 'bipartitionsBranchLabels' not in i: continue
		infile = open(inDIR+i, "r")
		ID = i.split('RAxML_bipartitionsBranchLabels.')[1]
		tre = infile.readline()
		boot = avg_bootstrap(tre)
		if boot >= cutoff:
			cmd = "cp "+inDIR+"RAxML_bestTree."+ID+" "+outDIR
			os.system(cmd)
		else:
			pass

		infile.close()




