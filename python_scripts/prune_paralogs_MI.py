"""
Input: homolog trees
Output: individual orthologs trees

if a tip is longer than the LONG_TIP_CUTOFF
and also long than 10 times its sister, cut it off
This is to fix the leftover trees that frequently has some long tips in it

If not to output 1-to-1 orthologs, for example, already analysed these
set OUTPUT_1to1_ORTHOLOGS to False
"""
import newick3,phylo3,os,sys
import trim_tips_module

REQUIRED = ['JA0456', 'JA0450', 'JA0816', 'JA0798', 'JA0723', 'JA0719', 'JA0711', 'JA0702', 'JA0701', 'JA0694', 'JA0432', 'JA0608', 'JA0726', 'Solyc']
OUTPUT_1to1_ORTHOLOGS = True 

def get_name(name):
	if '@' in name:
		return name[:6]
	else:
		return name[:5]

def get_clusterID(filename):
	return filename.split(".")[0]
	
def get_front_labels(node):
	leaves = node.leaves()
	return [i.label for i in leaves]

def get_back_labels(node,root):
	all_labels = get_front_labels(root)
	front_labels = get_front_labels(node)
	return set(all_labels) - set(front_labels)
	
def get_front_score(node):
	front_labels = get_front_labels(node)
	num_labels = len(front_labels)
	num_taxa = len(set([get_name(i) for i in front_labels]))
	if num_taxa == num_labels:
		return num_taxa
	return -1
	
def get_back_score(node,root):
	back_labels = get_back_labels(node,root)
	num_labels = len(back_labels)
	num_taxa = len(set([get_name(i) for i in back_labels]))
	if num_taxa == num_labels:
		return num_taxa
	return -1
	
def remove_kink(node,curroot):
	if node == curroot and curroot.nchildren == 2:
		#move the root away to an adjacent none-tip
		if curroot.children[0].istip: #the other child is not tip
			curroot = phylo3.reroot(curroot,curroot.children[1])
		else: curroot = phylo3.reroot(curroot,curroot.children[0])
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

def prune(score_tuple,node,root,pp_trees):
	if score_tuple[0] > score_tuple[1]: #prune front
		print "prune front"
		pp_trees.append(node)
		par = node.prune()
		if par != None and len(root.leaves()) >= 3:
			par,root = remove_kink(par,root)
		return root,node == root
	else:
		if node != root: #prune back
			par = node.parent #par--node<
			par.remove_child(node)
			if par.parent != None:
				par,root = remove_kink(par,root)
		node.prune()
		print "prune back"
		pp_trees.append(root)
		if len(node.leaves()) >= 3:
			node,newroot = remove_kink(node,node)
		else:
			newroot = node
		return newroot,False #original root was cutoff, not done yet
			

if __name__ == "__main__":
	if len(sys.argv) != 8:
		print "python prune_paralogs_MI.py homoTreeDIR tree_file_ending relative_tip_cutoff absolute_tip_cutoff1 absolute_tip_cutoff2 MIN_TAXA outDIR"
		print "LONG_TIP_CUTOFF is typically same value of the previous LONG_TIP_CUTOFF"
		sys.exit(0)

	inDIR = sys.argv[1]+"/"
	tree_file_ending = sys.argv[2]
	relative_tip_cutoff,absolute_tip_cutoff1,absolute_tip_cutoff2 = float(sys.argv[3]),float(sys.argv[4]),float(sys.argv[4])
	MIN_TAXA = int(sys.argv[6])
	outDIR = sys.argv[7]+"/"

	for i in os.listdir(inDIR):
		if not i.endswith(tree_file_ending): continue
		print i
		with open(inDIR+i,"r") as infile: #only 1 tree in each file
			intree = newick3.parse(infile.readline())
		curroot = intree
		pp_trees = []
		
		if get_front_score(curroot) >= MIN_TAXA: #No need to prune
			print "No pruning needed"
			if OUTPUT_1to1_ORTHOLOGS:
				os.system("cp "+inDIR+i+" "+outDIR+get_clusterID(i)+"_1to1ortho.tre")
		else: #scoring the tree
			going = True
			pp_trees = []
			while going: #python version of do..while loop
				highest = 0
				highest_node = None 
				score_hashes = {} #key is node, value is a tuple (front_score,back_score)
				for node in curroot.iternodes():
					front_score = get_front_score(node)
					back_score = get_back_score(node,curroot)
					score_hashes[node] = (front_score,back_score)
					if front_score > highest or back_score > highest:
						highest_node = node
						highest = max(front_score,back_score)
				if highest >= MIN_TAXA: #prune
					curroot,done = prune(score_hashes[highest_node],highest_node,curroot,pp_trees)
					if done or len(curroot.leaves()) < MIN_TAXA:
						going = False
						break
				else:
					going = False
					break
		
		if len(pp_trees) > 0:
			count = 1
			for tree in pp_trees:
				tree = trim_tips_module.trim(tree,relative_tip_cutoff,absolute_tip_cutoff1,absolute_tip_cutoff2)
				for x in REQUIRED:    # required species must be in the tree
					if x not in tree: continue
				if len(tree.leaves()) >= MIN_TAXA:
					with open(outDIR+get_clusterID(i)+"_MIortho"+str(count)+".tre","w") as outfile:	
						outfile.write(newick3.tostring(tree)+";\n")
					count += 1
