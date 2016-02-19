"""
Extract clades of sufficient ingroup coverages
If no outgroup, only take ingroup clades with no taxon repeats
If outgroup present, extract rooted ingroup clades and prune paralogs
"""

import phylo3,newick3,os,sys

HOMOTREE_ENDING = ".mm"
OUTGROUPS = ["Solyc"]
INGROUPS = ["JA0432","JA0694","JA0701","JA0014","JA0702","JA0711","JA0719","JA0723","JA0816","JA0798"]
MIN_INGROUP_TAXA = 10
#if pattern changes, change it here
#given tip label, return taxon name identifier
def get_name(label):
	if '@' in label:
		return label.replace("_R_","")[:6]
	else:
		return label.replace("_R_","")[:5]
		
def get_front_labels(node):
	leaves = node.leaves()
	return [i.label for i in leaves]

def get_back_labels(node,root):
	all_labels = get_front_labels(root)
	front_labels = get_front_labels(node)
	return set(all_labels) - set(front_labels) #labels do not repeat
	
def get_front_names(node): #may include duplicates
	labels = get_front_labels(node)
	return [get_name(i) for i in labels]

def get_back_names(node,root): #may include duplicates
	back_labels = get_back_labels(node,root)
	return [get_name(i) for i in back_labels]
	
def count_cary_names(node): #only count names with @
	labels = get_front_labels(node)
	cary_names = []
	for label in labels:
		if "@" in label:
			cary_names.append(get_name(label))
	return len(set(cary_names))
	
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
	
#input is a tree with both ingroups and more than 1 outgroups
def extract_ingroup_clades(root):
	print "extracting clades"
	inclades = []
	while True:
		max_score,direction,max_node = 0,"",None
		for node in root.iternodes():
			if node == root: continue
			front,back = 0,0
			front_names_set = set(get_front_names(node))
			for name in front_names_set:
				if name in OUTGROUPS:
					front = -1
					break
				else: front += 1
			back_names_set = set(get_back_names(node,root))
			for name in back_names_set:
				if name in OUTGROUPS:
					back = -1
					break
				else: back += 1
			if front > max_score:
				max_score,direction,max_node = front,"front",node
			if back > max_score:
				max_score,direction,max_node = back,"back",node
		#print max_score,direction
		if max_score >= MIN_INGROUP_TAXA:
			if direction == "front":
				inclades.append(max_node)
				kink = max_node.prune()
				if len(root.leaves()) > 3:
					newnode,root = remove_kink(kink,root)
				else: break
			elif direction == "back":
				par = max_node.parent
				par.remove_child(max_node)
				max_node.prune()
				inclades.append(phylo3.reroot(root,par))#flip dirction
				if len(max_node.leaves()) > 3:
					max_node,root = remove_kink(max_node,max_node)
				else: break
		else: break
	return inclades
	
def get_ortho_from_rooted_inclade(inclade):
	orthologs = [] #store ortho clades
	clades = [inclade]
	while True:
		newclades = [] #keep track of subclades generated in this round
		for clade in clades:
			num_taxa = len(set(get_front_names(clade)))
			num_tips = len((get_front_labels(clade)))
			if count_cary_names(clade) < MIN_INGROUP_TAXA:
				pass #not enough taxa
			elif num_taxa == num_tips and count_cary_names(clade) >= MIN_INGROUP_TAXA:
				orthologs.append(clade) #enough taxa and all taxa are unique
			else: #enough taxa but has duplicated taxa
				for node in clade.iternodes(order=0): #PREORDER, root to tip
					if node.istip: continue
					#traverse the tree from root to tip
					child0,child1 = node.children[0],node.children[1]
					name_set0 = set(get_front_names(child0))
					name_set1 = set(get_front_names(child1))
					if len(name_set0.intersection(name_set1)) > 0:
						if node == clade:
							newclades += [child0,child1] #break by bifid at the base
						elif len(name_set0) > len(name_set1): #cut the side with less taxa
							node.remove_child(child1)
							child1.prune()
							node,clade = remove_kink(node,clade) #no rerooting here
							newclades += [clade,child1]
						else:
							node.remove_child(child0)
							child0.prune()
							node,clade = remove_kink(node,clade) #no rerooting here
							newclades += [clade,child0]
						break
		if newclades == []: break
		clades = newclades
	return orthologs


if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "python prune_paralogs_RT.py homoTreeDIR outDIR"
		sys.exit(0)
	
	homoDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"

	for treefile in os.listdir(homoDIR):
		if treefile[-len(HOMOTREE_ENDING):] != HOMOTREE_ENDING: continue
		print treefile
		with open(homoDIR+treefile,"r") as infile:
			 intree = newick3.parse(infile.readline())
		curroot = intree
		
		#check whehter there's any outgroup at all
		all_names = get_front_names(curroot)
		outgroup_presnet = False
		for name in all_names:
			if name not in INGROUPS and name not in OUTGROUPS:
				print "check name",name
				sys.exit()
			if name in OUTGROUPS:
				outgroup_presnet = True
				break
		
		if outgroup_presnet: #at least one outgroup present, root and cut inclades
			inclades = extract_ingroup_clades(curroot)
			print len(inclades),"inclades extracted"
			inclade_count = 0
			for clade in inclades:
				inclade_count += 1
				inclade_name = outDIR+treefile.replace(HOMOTREE_ENDING,".inclade")+str(inclade_count)
				with open(inclade_name,"w") as outfile:
					outfile.write(newick3.tostring(clade)+";\n")
				orthologs = get_ortho_from_rooted_inclade(clade)
				ortho_count = 0
				for ortho in orthologs:
					ortho_count += 1
					with open(inclade_name+".ortho"+str(ortho_count)+".tre","w") as outfile:
						outfile.write(newick3.tostring(ortho)+";\n")
		else: #no outgroup
			#only output tree when there is no taxon repeats
			#do not attempt to infer direction of gene duplication without outgroup info
			num_taxa = len(all_names)
			if len(all_names) == len(set(all_names)):
				if count_cary_names >= MIN_INGROUP_TAXA:
					with open(outDIR+treefile.replace(HOMOTREE_ENDING,".unrooted-inclade"),"w") as outfile:
						outfile.write(newick3.tostring(curroot)+";\n")
					with open(outDIR+treefile.replace(HOMOTREE_ENDING,".unrooted-ortho.tre"),"w") as outfile:
						outfile.write(newick3.tostring(curroot)+";\n")
				else: print "only",num_taxa,"taxa present in unrooted tree"
			else: print "duplicated taxa in unrooted tree"
		continue
		

