import sys,os,newick3,phylo3

"""
start from fasttree results from the first round of alignment

not automating the first round so that one can check the alignments and trees
for problematic clusters and process clusters that are unusually large
and adjust filters for all-by-all blast

Walk through the tree and cut branches longer than the first cutoff, re-align until
no branch is longer than the first cutoff

then cut branches longer than the second cutoff
re-align until no branch is longer than the second cutoff

write new fasta files for prank
"""

MIN_INGROUP_TAXA = 10
#all ingroup datasets have "@" after the four-character taxon name code
#will skip fasta file if the number of uniq taxon codes are less than this

NUM_MAFFT_CORES = 8
#number of processors mafft uses
#generally should be less than 10

cutoff1 = 0.3
cutoff2 = 0.1

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

def count_ingroups(node):
	labels = get_leaf_labels(node.leaves())
	names = []
	for label in labels:
		name = get_name(label)
		if "@" in label and name not in names:
			names.append(name)
	return len(names)

#do not score terminal branches
#long terminal branches usually do not disturb the alignment as much
def find_longest_internal_branch_length(tree):
	longest_internal_branch_length = 0.0
	for node in tree.iternodes():
		if node != tree and not node.istip:
			longest_internal_branch_length = max(longest_internal_branch_length,node.length)
	return longest_internal_branch_length

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
def cut_long_branches(curroot,cutoff):
	going = True
	subtrees = [] #store all subtrees after cutting
	if curroot.nchildren == 2: #fix the root
		#move the root away to an adjacent none-tip internal node
		if curroot.children[0].istip: #the other child is not tip
			curroot = phylo3.reroot(curroot,curroot.children[1])
		else: #tree has >=4 leaves so the other node cannot be tip
			curroot = phylo3.reroot(curroot,curroot.children[0])
	while going:
		going = False #only keep going if long branches were found during last round
		for node in curroot.iternodes(): #Walk through nodes
			if node != curroot and node.length > cutoff:
				subtrees.append(node)
				node = node.prune()
				if len(curroot.leaves()) >= 4:
					node,curroot = remove_kink(node,curroot)
					going = True
				break
	subtrees.append(curroot) #write out the residue after cutting
	return subtrees

def mafft_align(fasta_file):
	count = 0 #record how many sequences in the fasta file
	with open(fasta_file,"r") as infile:
		for line in infile:
			if line[0] == ">": count += 1
	if count >= 1000: alg = "--auto" #so that the run actually finishes!
	else: alg = "--genafpair --maxiterate 1000"
	com = "mafft "+alg+" --nuc --thread 8"
	com += str(NUM_MAFFT_CORES)+" "+fasta_file+" > "+fasta_file+".aln"
	print com
	os.system(com)
	return fasta_file+".aln"

#give the path of the alignment file
#remove columns with occupancy lower than 0.1
#remove seqs shorter than 200 bp after filter columns
def phyutility_clean_alignment(alignment):
	cmd = "./phyutility -clean 0.1 -in "+alignment+" -out "+alignment+"-pht"
	print cmd
	os.system(cmd)	
	seqid = ""
	first = True
	infile = open(alignment+"-pht","r")
	outfile = open(alignment+"-cln","w")
	for line in infile:
		line = line.strip()
		if len(line) == 0: continue #skip empty lines
		if line[0] == ">":
			if seqid != "": #not at the first seqid
				if len(seq.replace("-","")) >= 200:
					outfile.write(">"+seqid+"\n"+seq+"\n")
			seqid,seq = line.split(" ")[0][1:],""
		else: seq += line.strip()
	#process the last seq
	if len(seq.replace("-","")) >= 200:
		outfile.write(">"+seqid+"\n"+seq+"\n")
	infile.close()
	outfile.close()
	return alignment+"-cln"

def fasttree(cleaned_alignment):
	com = "FastTreeMP -wag "+cleaned_alignment+" >"+cleaned_alignment+".tre"
	print com
	os.system(com)
	return cleaned_alignment+".tre"

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "python cut_long_branches_iter.py inDIR outDIR >log"
		sys.exit(0)
	
	DIR = sys.argv[1]+"/"
	for i in os.listdir(DIR): #go through fasta files in the input directory
		if i[-9:] != ".fasttree": continue
		fasta_file = DIR+i.replace(".aln-cln.fasttree","")
		tree_file = DIR+i
		with open(tree_file,"r") as infile:
			first_line = infile.readline() #there's only 1 tree in each file
		if first_line.strip() == "": continue #empty file after trimming
		intree = newick3.parse(first_line.replace ("-","_"))
		if count_ingroups(intree) < MIN_INGROUP_TAXA:
			continue #skip trees with few ingroups
		ccID = i.split(".")[0] #looks like cc9
		trees = [intree]
		
		#if intree has no long branches at all jus use the original fasta and alignment
		longest_branch_length = 0.0
		for node in intree.iternodes():
			if node != intree:
				longest_branch_length = max(longest_branch_length,node.length)
		if longest_branch_length < cutoff2:
			os.system("cp "+fasta_file+" "+fasta_file.replace(".fa",".to_prank.fa"))
			os.system("cp "+fasta_file+".aln-cln "+fasta_file.replace(".fa",".to_prank.fa.aln-cln"))
			continue

		#read fasta file
		infile = open(fasta_file,"r")
		seqDICT = {} #key is seq id, value is the sequence
		seqid = ""
		for line in infile:
			line = line.strip()
			print line
			if len(line) > 0:
				if line[0] == ">":
					if seqid != "":
						seqDICT[seqid] = seq
					seqid,seq = line[1:].replace("-","_"),""
				else: seq += line.strip()
		seqDICT[seqid] = seq #add the last record
		infile.close()

		#cut by the first cutoff
		print i,"Cutting branches longer than",cutoff1
		outfile = open(DIR+ccID+".cut1.trees","w")
		newtrees = [] #store trees in need of cutting for the next round
		count = 0
		while True:
			for tree in trees:
				if find_longest_internal_branch_length(tree) < cutoff1:
					outfile.write(newick3.tostring(tree) +";\n")
					#can be the original tree or cut tree. no need to cut tip here
				else:
					subtrees = cut_long_branches(tree,cutoff1)
					for subtree in subtrees:
						if count_ingroups(subtree) < MIN_INGROUP_TAXA: continue
						count += 1
						newname = DIR+ccID+".cut1-"+str(count)
						with open(newname+".cutbranch","w") as outfile1: #record the cut branch
							outfile1.write(newick3.tostring(subtree)+";\n")
						with open(newname+".fa","w") as outfile2: #output fasta
							for label in get_leaf_labels(subtree.leaves()):
								outfile2.write(">"+label+"\n"+seqDICT[label]+"\n")
						newaln = mafft_align(newname+".fa")
						newcln = phyutility_clean_alignment(newaln)
						newtreefile = fasttree(newcln)
						with open(newtreefile) as infile: #update the tree
							newtrees.append(newick3.parse(infile.readline()))
			if newtrees == []: break
			trees = newtrees
			newtrees = []
		outfile.close()

		trees = [] #update the tree record for the next round
		with open(DIR+ccID+".cut1.trees","r") as infile:
			for line in infile:
				if len(line) > 3:
					trees.append(newick3.parse(line))
		
		#cut by the second cutoff
		print len(trees),"trees cutting by",cutoff2
		outfile = open(DIR+ccID+".cut2-final.trees","w")
		count = 0
		newtrees = [] #store trees in need of cutting for the next round
		while True:
			for tree in trees:
				print find_longest_internal_branch_length(tree)
				if find_longest_internal_branch_length(tree) < cutoff2:
					subtrees = cut_long_branches(tree,cutoff2)
					for subtree in subtrees:
						if count_ingroups(subtree) >= MIN_INGROUP_TAXA:
							outfile.write(newick3.tostring(subtree)+";\n")
				else:
					subtrees = cut_long_branches(tree,cutoff2)
					print len(subtrees)
					for subtree in subtrees:
						if count_ingroups(subtree) < MIN_INGROUP_TAXA: continue
						count += 1
						newname = DIR+ccID+".cut2-"+str(count)
						with open(newname+".cutbranch","w") as outfile1:
							outfile1.write(newick3.tostring(subtree)+";\n")
						with open(newname+".fa","w") as outfile2: #output fasta
							for label in get_leaf_labels(subtree.leaves()):
								outfile2.write(">"+label+"\n"+seqDICT[label]+"\n")
						newaln = mafft_align(newname+".fa")
						newcln = phyutility_clean_alignment(newaln)
						newtreefile = fasttree(newcln)
						with open(newtreefile) as infile: #update the tree
							newtrees.append(newick3.parse(infile.readline()))
			if newtrees == []: break
			trees = newtrees
			newtrees = []
		outfile.close()
		
		try:
			with open(DIR+ccID+".cut2-final.trees","r") as infile:
				count = 1
				for line in infile:
					if len(line) <= 3: continue
					tree = newick3.parse(line)
					newfasta = sys.argv[2]+"/"+ccID+"-"+str(count)+".final.fa"
					with open(newfasta,"w") as outfile2: #output fasta
						for label in get_leaf_labels(tree.leaves()):
							outfile2.write(">"+label+"\n"+seqDICT[label]+"\n")
					count += 1
		except: pass

