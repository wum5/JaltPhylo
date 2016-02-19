import sys,os
import phylo3,newick3
from tree_utils import get_name

"""
to change sequence names with taxa names and make trees more readable,
and output a new file named infile.names

Create a tabular file that each line contains
code	taxon_name
separated by tab
"""

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "python taxon_name_subst.py table treefile"
		sys.exit(0)
	
	DICT = {} #key is seq acronym, value is full taxon name, separated by tab
	with open(sys.argv[1], "rU") as infile:
		for line in infile:
			spls = line.strip().split("\t")
			if len(spls) > 1:
				DICT[spls[0]] = spls[1]
	print DICT
	
	treefile = sys.argv[2]
	with open(treefile,"r") as infile:
		intree = newick3.parse(infile.readline())
	for i in intree.iternodes():
		if i.istip:
			taxonID = get_name(str(i.label))
			i.label = taxonID+"_"+DICT[taxonID]

	with open(treefile+".name","w") as outfile:
		outfile.write(newick3.tostring(intree)+";\n")