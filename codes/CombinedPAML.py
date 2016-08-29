import sys

ID1, ID2, ID3 = [], [], []
function, NON, SYN = [], [], []
NSratio1, NSratio2, ltr = [], [], []
size, taxon = [], []


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "usage: python CombinedPAML.py NSFile PAMLFile FuncFile"
		sys.exit()

	NSFile = open(sys.argv[1],"r") 
	for line in NSFile:
		if 'Solyc' not in line: continue
		line = line.rstrip()
		ID2.append(line.split('\t')[0])
		NON.append(line.split('\t')[1])
		SYN.append(line.split('\t')[2])
		NSratio1.append(line.split('\t')[3])


	PAMLFile = open(sys.argv[2],"r")
	for line in PAMLFile:
		if 'Solyc' not in line: continue
		line = line.rstrip()
		ID3.append(line.split('\t')[0])
		taxon.append(line.split('\t')[1])
		size.append(line.split('\t')[2])
		ltr.append(line.split('\t')[3])
		NSratio2.append(line.split('\t')[22])


	FuncFile = open(sys.argv[3],"r") 
	for line in FuncFile:
		if 'Solyc' not in line: continue
		line = line.rstrip()
		ID1.append(line.split('\t')[0])
		function.append(line.split('\t')[1])


	print "Genes\tSize\tNON\tSYN\tNON/SYN\tltrscore\tPAML dN/dS\tAnnotations"
	for x in xrange(len(ID2)):
		currFunc = 'Unknown function'
		for y in xrange(len(ID3)):
			if ID2[x] == ID3[y]:
				for z in xrange(len(ID1)):
					if ID1[z] in ID2[x] :
						currFunc = function[z]
						break
					else:
						''
				break

		print ID2[x]+'\t'+size[y]+'\t'+NON[x]+'\t'+SYN[x]+'\t'+NSratio1[x]+'\t'+ltr[y]+'\t'+NSratio2[y]+'\t'+currFunc 


