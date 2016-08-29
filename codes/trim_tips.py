"""
Trim tips that sticking out (> relative_cutoff and >10 times longer than sister; here just for ingroup species)
Also trim any tips that are > absolute_cutoff1 (outgroup) and  absolute_cutoff2 (ingroup)
"""

import newick3,phylo3,os,sys
from Bio import SeqIO

#if taxon id pattern changes, change it here
def get_name(name):
	if '@' in name:
		return name[:6]
	else:
		return name[:5]

	
#return the outlier tip, with abnormal high contrast and long branch
def check_countrast_outlier(node0,node1,above0,above1,id0,id1):
	global relative_cutoff
	if 'Solyc' in str(id0) or 'Solyc' in str(id1):
		return None
	if node0.istip and above0>relative_cutoff:
		if above1 == 0.0 or above0/above1>10:
			return node0
	if node1.istip and above1>relative_cutoff:
		if above0 == 0.0 or above1/above0>10 :
			return node1
	return None
	
def cut_long_tips(curroot):
	to_remove = []
	global absolute_cutoff1, absolute_cutoff2
	for i in curroot.iternodes(order=1):#POSTORDER
		if i.nchildren == 0: #at the tip
			i.data['len'] = i.length
			if i.length > absolute_cutoff1 and 'Solyc' in str(i.label):
				to_remove.append(i)
			if i.length > absolute_cutoff2 and 'JA0' in str(i.label):
				to_remove.append(i)
		elif i.nchildren == 1:
			print "kink in tree"
			sys.exit()
		elif i.nchildren == 2: #bifurcating internal nodes
			child0,child1 = i.children[0],i.children[1]
			above0,above1 = child0.data['len'],child1.data['len']
			id0, id1 = child0.label, child1.label
			i.data['len'] = ((above0+above1)/2.)+i.length #stepwise average
			outlier = check_countrast_outlier(child0,child1,above0,above1,id0,id1)
			if outlier != None:
				to_remove.append(outlier)
		else: #3 or more branches from this node
			total_len = 0
			nchild = i.nchildren
			for child in i.children:
				total_len += child.data['len']
			i.data['len'] = total_len / float(i.nchildren)
			#print "total_len:",total_len,"nchildren:",i.nchildren,"ave:",i.data['len']
			for index1 in range(nchild): #do all the pairwise comparison
				for index2 in range(nchild):
					if index2 <= index1: continue
					child1,child2 = i.children[index1],i.children[index2]
					id1, id2 = child1.label, child2.label
					above1,above2 = child1.data['len'], child2.data['len']
					outlier = check_countrast_outlier(child1,child2,above1,above2,id1,id2)
					if outlier != None:
						to_remove.append(outlier)
	to_remove = list(set(to_remove))


	if len(to_remove) > 0:
		for node in to_remove:
			print node.label,node.length
			node = node.prune()
			if len(curroot.leaves()) > 3:
				node,curroot = remove_kink(node,curroot)
			else:
				print "Less than four tips left"
	return curroot
	
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

	
if __name__ == "__main__":
	if len(sys.argv) != 7:
		print "python trim_tips.py treDIR tree_file_ending outDIR relative_cutoff absolute_cutoff1 absolute_cutoff2"
		sys.exit(0)

	treDIR = sys.argv[1]+"/"
	file_ending = sys.argv[2]
	outDIR = sys.argv[3]+"/"
	relative_cutoff = float(sys.argv[4])
	absolute_cutoff1 = float(sys.argv[5])
	absolute_cutoff2 = float(sys.argv[6])
	
	done = [] #record clusterIDs that are done
	for i in os.listdir(treDIR):
		if i[-3:] == ".tt":
			done.append(i.split(".")[0])
	print done
	
	filecount = 0
	l = len(file_ending)
	for i in os.listdir(treDIR):
		if file_ending in i:
			clusterID = i.split(".")[0]
			if clusterID in done: continue
			print i
			filecount += 1
			with open(treDIR+i,"r") as infile:
				intree = newick3.parse(infile.readline())
			with open(outDIR+i+".tt","w") as outfile:
				outfile.write(newick3.tostring(cut_long_tips(intree))+";\n")
	
	if filecount == 0:
		print "No file name with",file_ending,"found in the treDIR"