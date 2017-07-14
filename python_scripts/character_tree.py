from Bio import Phylo
from cStringIO import StringIO
import sys, os


# lineages with derived traits
#consensus = set(['J.quipuscoae','J.calliantha','J.yungayensis','J.biflora','J.sinuosa','J.aijana','J.umbellata','J.grandibaccata','J.dendroidea','J.incahuasina'])
#nectar = set(['J.quipuscoae','J.calliantha','J.biflora','J.aijana','J.umbellata','J.grandibaccata','J.dendroidea','J.incahuasina'])
campan = set(['JA0798','JA0711','JA0010','JA0719'])
campan1 = set(['JA0798','JA0711'])
campan2 = set(['JA0010','JA0719'])


def ancestral_brcLen(substr):
	n = 1
	for x in range(len(substr)):
		if substr[x] == '(':
			n += 1
		elif substr[x] == ')':
			n -= 1
		if n == 0:
			return substr[x+2:x+24]

if __name__ == "__main__":
	
	if len(sys.argv) != 3:
		print "usage: python character_tree.py inDIR IDfile"
		sys.exit()
	
	inDIR = sys.argv[1]+'/'		
	IDfile = open(sys.argv[2], "r")	
	geneID, geneFct = [], []
	outfile1 = open("single_origin.txt", "w")
	outfile2 = open("double_origin.txt", "w")
	
	for line in IDfile:
		line = line.rstrip()
		geneID.append(line.split('\t')[0])
		geneFct.append(line.split('\t')[1])
		
	for i in os.listdir(inDIR):
		if 'best' not in i: continue
		mynode, mynode1, mynode2 = False, False, False
		infile = open(inDIR+i, "r")
		line = infile.readline().rstrip()
		handle = StringIO(line)
		tree = Phylo.read(handle, "newick")

		for node in tree.get_nonterminals():
			subnodes = set([x.name for x in node.get_terminals()])
			if subnodes == campan1:
				mynode1 = True
			elif subnodes == campan2:
				mynode2 = True
			elif subnodes == campan:
				mynode = True
		
		if mynode == True:
			for x in xrange(len(geneID)):
				if geneID[x] in i:
					outfile1.write("%s\t%s\t%s\n" % (geneID[x],geneFct[x],line))
					break
				
		elif mynode1 == True and mynode2 == True and mynode == False:
			start = line.find('JA0711')
			if line[start+29] == ')':
				BrcLen1 = line[start+31:start+53]
				#print BrcLen1
				BrcLenA1 = ancestral_brcLen(line[start+53:])
				#print BrchLenA1			
			else:
				BrcLen1 = line[start+61:start+83]
				#print BrcLen1
				BrcLenA1 = ancestral_brcLen(line[start+83:])
				#print BrchLenA1
						
			start = line.find('JA0719')
			if line[start+29] == ')':
				BrcLen2 = line[start+31:start+53]
				#print BrcLen2
				BrcLenA2 = ancestral_brcLen(line[start+53:])
				#print BrchLenA2				
			else:
				BrcLen2 = line[start+61:start+83]
				#print BrcLen2
				BrcLenA2 = ancestral_brcLen(line[start+83:])
				#print BrchLenA2

			print BrcLenA1, BrcLenA2
			if float(BrcLen1) < 0.00001 or float(BrcLen2) < 0.00001: continue
			if float(BrcLenA1) < 0.00001 or float(BrcLenA2) < 0.00001: continue
			
			#print line
			for x in xrange(len(geneID)):
				if geneID[x] in i:
					outfile2.write("%s\t%s\t%s\n" % (geneID[x],geneFct[x],line))
					break
		
		infile.close()	
		
	outfile1.close()
	outfile2.close()
	
	