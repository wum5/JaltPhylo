"""
Separate subtrees connected by long branches
Input trees: tree files from tip masking
"""

import newick3,phylo3,os,sys,math
from Bio import SeqIO

#if taxon id pattern changes, change it here
def get_name(name):
	if '@' in name:
		return name[:6]
	else:
		return name[:5]
	
def get_leaf_labels(leaves):
	labels = []
	for i in leaves:
		labels.append(i.label)
	return labels

def count_ingroup_taxa(node):
	labels = get_leaf_labels(node.leaves())
	names = []
	for label in labels:
		if "@" in label:
			names.append(get_name(label))
	return len(set(names))
	
def count_outgroup_taxa(node):
	labels = get_leaf_labels(node.leaves())
	names = []
	for label in labels:
		if "@" not in label:
			names.append(get_name(label))
	return len(set(names))
	
#smooth the kink created by prunning
#to prevent creating orphaned tips after prunning twice at the same node
def remove_kink(node,curroot):
	if node == curroot and curroot.nchildren == 2:
		#move the root away to an adjacent none-tip internal node
		if curroot.children[0].istip: #the other child is not tip
			curroot = phylo3.reroot(curroot,curroot.children[1])
		else: #tree has >=4 leaves so the other node cannot be tip
			curroot = phylo3.reroot(curroot,curroot.children[0])
	#---node---< all nodes should have one child only now
	length = node.length + (node.children[0]).length
	par = node.parent
	kink = node
	node = node.children[0]
	#parent--kink---node<
	par.remove_child(kink)
	par.add_child(node)
	node.length = length
	return node,curroot

#cut long branches and output all subtrees regardless of size
def cut_long_internal_branches(curroot,cutoff):
	going = True
	subtrees = [] #store all subtrees after cutting
	while going:
		going = False #only keep going if long branches were found during last round
		for node in curroot.iternodes(): #Walk through nodes
			if node.istip or node == curroot: continue
			child0,child1 = node.children[0],node.children[1]
			if node.length > cutoff:
				print node.length
				if not child0.istip and not child1.istip and child0.length+child1.length>cutoff:
					print child0.length + child1.length
					if count_ingroup_taxa(child0) >= 4:
						subtrees.append(child0)
					if count_ingroup_taxa(child1) >= 4:
						subtrees.append(child1)						
				else: subtrees.append(node)
				node = node.prune()
				if len(curroot.leaves()) > 2: #no kink if only two left
					node,curroot = remove_kink(node,curroot)
					going = True
				break
	if count_ingroup_taxa(curroot) >= min_ingroup_taxa:
		subtrees.append(curroot) #write out the residue after cutting
	return subtrees
	

if __name__ == "__main__":
	if len(sys.argv) != 6:
		print "python cut_long_internal_branches.py inDIR internal_branch_length_cutoff minimal_ingroup_taxa minimal_outgroup_taxa outDIR"
		sys.exit(0)

	inDIR = sys.argv[1]+"/"
	branch_len_cutoff = float(sys.argv[2]) #cut branches longer than this
	min_ingroup_taxa = int(sys.argv[3])
	min_outgroup_taxa = int(sys.argv[4])
	outDIR = sys.argv[5]+"/"
	
	filecount = 0
	for i in os.listdir(inDIR):
		if "best" not in i and i[-3:] != ".mm": continue
		filecount += 1
		if "best" in i:
			clusterID = i.split(".")[1]
		else: clusterID = i.split(".")[0]
		print clusterID
		with open(inDIR+i,"r") as infile: #only 1 tree in each file
			intree = newick3.parse(infile.readline())
		curroot = intree
		if count_ingroup_taxa(curroot) < min_ingroup_taxa: continue
		subtrees = cut_long_internal_branches(curroot,branch_len_cutoff)
		if len(subtrees) > 0:
			count = 1
			for subtree in subtrees:
				if count_ingroup_taxa(subtree)>=min_ingroup_taxa and count_outgroup_taxa(subtree)>=min_outgroup_taxa:
					#fix bifurcating roots from cutting
					if subtree.nchildren == 2:
						subtree,subtree = remove_kink(subtree,subtree)
					with open(outDIR+clusterID+"_subtree"+str(count)+".tre","w") as outfile:
						outfile.write(newick3.tostring(subtree)+";\n")
						count += 1
	
	if filecount == 0:
		print "No file end with",file_ending,"found"
			
		
		
