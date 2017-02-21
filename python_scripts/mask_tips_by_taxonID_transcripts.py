"""
Input is a dir of trees

Mask both mono- and paraphyletic tips that belong to the same taxon
Keep the tip that has the most un-ambiguous, well-aligned charactors in the trimmed alignment

Change ALIGNMENT_FILE_ENDING to match the corresponding alignment files
"""

import newick3,phylo3,os,sys
from Bio import SeqIO

#ALIGNMENT_FILE_ENDING = ".sate.aln-cln"
ALIGNMENT_FILE_ENDING = ".final.fa.aln-cln"
#ALIGNMENT_FILE_ENDING = ".aln-cln"
#if taxon id pattern changes, change it here
def get_name(name):
	if '@' in name:
		return name[:6]
	else:
		return name[:5]

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

def get_sister_tips(tip):
	sisters = tip.get_sisters()
	sister_tips = []
	for sister in sisters:
		if sister.istip: sister_tips.append(sister)
	print tip.label, len(sister_tips)
	return sister_tips

def monophyly_masking(curroot,unamb_chrDICT):
	going = True
	while going and len(curroot.leaves()) >= 4:
		going = False
		for node in curroot.iternodes(): #walk through nodes
			if not node.istip: continue #only look at tips
			for sister in node.get_sisters():
				if sister.istip and get_name(node.label)==get_name(sister.label): #masking
					#print node.label,unamb_chrDICT[node.label],sister.label,unamb_chrDICT[sister.label]
					if unamb_chrDICT[node.label] > unamb_chrDICT[sister.label]:
						node = sister.prune()			
					else: node = node.prune()
					if len(curroot.leaves()) >= 4:
						if (node==curroot and node.nchildren==2) or (node!=curroot and node.nchildren==1):
							node,curroot = remove_kink(node,curroot)
					going = True
					break
	return curroot
	
def paraphyly_masking(curroot,unamb_chrDICT):
	going = True
	while going and len(curroot.leaves()) >= 4:
		going = False
		for node in curroot.iternodes(): #walk through nodes
			if not node.istip: continue #only look at tips
			parent = node.parent
			if node == curroot or parent == curroot:
				continue #no paraphyletic tips for the root
			for para in parent.get_sisters():
				if para.istip and get_name(node.label)==get_name(para.label):
					if unamb_chrDICT[node.label] > unamb_chrDICT[para.label]:
						node = para.prune()
					else: node = node.prune()
					if len(curroot.leaves()) >= 4:
						if (node==curroot and node.nchildren==2) or (node!=curroot and node.nchildren==1):
							node,curroot = remove_kink(node,curroot)
					going = True
					break
	return curroot
	
if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "python mask_tips_by_taxonID_transcripts.py treDIR aln-clnDIR outDIR mask_para(y/n)"
		sys.exit(0)

	treDIR = sys.argv[1]+"/"
	clnDIR = sys.argv[2]+"/"
	outDIR = sys.argv[3]+"/"
	if sys.argv[4] == "y": mask_para = True
	elif sys.argv[4] == "n": mask_para = False
	else:
		print "mask_para? y/n"
		sys.exit()
	filecount = 0
	for i in os.listdir(treDIR):
		if i[-3:] == ".tt" and i[-3:] != ".mm":
			with open(treDIR+i,"r") as infile:
				intree = newick3.parse(infile.readline())
			print i
			clusterID = i.split("_")[0]
			filecount += 1
			unamb_chrDICT = {} #key is seqid, value is number of unambiguous chrs
			with open(clnDIR+clusterID+ALIGNMENT_FILE_ENDING) as handle:
				for record in SeqIO.parse(handle,"fasta"):
					seqid,seq = str(record.id), str(record.seq)
					for ch in ['-','X',"x","?","*"]:
						seq = seq.replace(ch,"") #ignore gaps, xs and Xs
					unamb_chrDICT[seqid] = len(seq)
			curroot = monophyly_masking(intree,unamb_chrDICT)
			if mask_para:
				curroot = paraphyly_masking(curroot,unamb_chrDICT)
			with open(outDIR+i+".mm","w") as outfile:
				outfile.write(newick3.tostring(curroot)+";\n")
	if filecount == 0:
		print "No file name with 'best' or 'tt' or 'fasttree' found in the treDIR"
