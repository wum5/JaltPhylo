"""
Filter the aligned and cleaned ortholog matrices by number of taxa and characters
Write out name of matrices that passed the filter
Also write out supermatrix stats
"""
from Bio import SeqIO
from tree_utils import get_name
import sys,os

MATRIX_FILE_ENDING = ".aln-cln"
MATRIX_FILE_HEADING = "cc"

if __name__ == "__main__":
	if len(sys.argv) != 6:
		print "usage: python concatenate_matrices.py aln-clnDIR numofsitesFilter numoftaxaFilter dna/aa outfile"
		sys.exit()

	clnDIR = sys.argv[1]+"/"
	sites_filter = int(sys.argv[2])
	taxa_filter = int(sys.argv[3])
	if sys.argv[4] == "aa":
		seqtype = " -aa"
		model = "WAG" 
	elif sys.argv[4] == "dna":
		seqtype = "" 	
		model = "DNA"
	else:
		print "Input data type: DNA or aa"
		sys.exit()
	outname = sys.argv[5]
	
	print "Filtering ortholog matrixes and writing the list of selected matrixes"
	outfile = open(outname+"_selected_ortho_matrixes","w")
	selected_matrixes = []
	for i in os.listdir(clnDIR):
		if i[-len(MATRIX_FILE_ENDING):] != MATRIX_FILE_ENDING: continue
		with open(clnDIR+i,"rU") as handle:
			seqs = list(SeqIO.parse(handle,"fasta"))
			num_taxa = len(seqs)
			if num_taxa >= taxa_filter:
				num_sites = len(seqs[0].seq)
				if num_sites >= sites_filter:
					outfile.write(str(i)+"\n")
					selected_matrixes.append(i)
	outfile.close()
	print len(selected_matrixes),"matrices passed the filter"
	print "List of selected ortholog matrixes written to",outname+"_selected_ortho_matrixes"
	
	print "Getting matrix occupancy stats"
	taxon_occupancy = {}
	#key is taxon name, value is [times present in a matrix,total length for this taxon]
	total_aligned_len = 0 #record how long the final concatenated matrix is
	cmd = "phyutility -concat"+seqtype+" -out "+outname+".nex -in "
	for i in selected_matrixes:
		cmd += clnDIR+i+" "
		handle = open(clnDIR+i,"rU")
		first = True
		for seq_record in SeqIO.parse(handle,"fasta"):
			taxon,seq = get_name(str(seq_record.id)),str(seq_record.seq)
			if first:
				total_aligned_len += len(seq)
				first = False
			if taxon not in taxon_occupancy:
				taxon_occupancy[taxon] = [0,0]
			taxon_occupancy[taxon][0] += 1
			taxon_occupancy[taxon][1] += len((seq.replace("-","")).replace("?",""))
		handle.close()
	cmd += "\n"
	
	total_ortho = len(selected_matrixes)
	outfile = open(outname+"_taxon_occupancy_stats","w")
	outfile.write("taxon\t#orthologs\t#total_charactors\tperc_orthologs\tperc_charactors\n")
	sum_char = 0
	for taxon in taxon_occupancy:
		times,chars = taxon_occupancy[taxon][0],taxon_occupancy[taxon][1]
		sum_char += chars
		out = taxon+"\t"+str(times)+"\t"+str(chars)+"\t"
		out += str(times/float(total_ortho))+"\t"+str(chars/float(total_aligned_len))+"\n"
		outfile.write(out)
	total_taxa = len(taxon_occupancy)
	out = "\nSupermatrix dimension "+str(total_taxa)+" taxa, "
	out += str(total_ortho)+" loci and "+str(total_aligned_len)+" aligned columns\n"
	out += "Overall matrix occupancy "+str(sum_char/float(total_taxa*total_aligned_len))+"\n"
	outfile.write(out)
	outfile.close()
	print "Supermatrix taxon occupancy stats written to",outname+"_taxon_occupancy_stats"
	print "Waiting for concatenation to finish. This may take several minutes..."
	with open(outname+".temp.sh","w") as outfile:
		outfile.write(cmd)
	os.system("sh "+outname+".temp.sh")
	os.system("rm "+outname+".temp.sh")
	
	#convert the .nex file to .phy and .model files for raxml
	infile = open(outname+".nex","r")
	outfile = open(outname+".phy","w")
	for line in infile:
		line = line.strip()
		if len(line) < 10: continue
		if line[0]=="#" or line[:5]=="BEGIN" or line[:6]=="MATRIX" or line=="END;" or line[:6]=="FORMAT":
			continue
		if line[0] == "[":
			line = line.replace("[",model+",")
			line = line.replace(" ]","")
			line = line.replace(" "+MATRIX_FILE_HEADING,"\n"+model+","+MATRIX_FILE_HEADING)
			line = line.replace(" ","=")
			with open(outname+".model","w") as outfile2:
				outfile2.write(line.strip()+"\n")
				#make sure that wc -l will get how many partitions
		elif line[:10]=="DIMENSIONS":
			ntax = (line.split("NTAX=")[1]).split(" ")[0]
			nchar = (line.split("NCHAR=")[1]).replace(";","")
			outfile.write(ntax+" "+nchar+"\n")
		else:
			spls=line.split("\t")
			outfile.write(spls[0]+" "+spls[1]+"\n")
	infile.close()
	outfile.close()
	print "outfiles written",outname+".phy",outname+".model"
	os.system("rm "+outname+".nex") #remove intermediate .nex file