import numpy, sys, os, time
from Bio import SeqIO


File_Ending = '.fa'

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print "usage: python lengthFilter.py inDIR outDIR"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	outDIR = sys.argv[2]+"/"

	start_time = time.time()	## record the start time
	
	for i in os.listdir(inDIR):
		if i[-len(File_Ending):] != File_Ending: continue
		handle = open(inDIR+i, "rU")
		infile = SeqIO.parse(handle, "fasta")	
		
		check = []
		unmasked_sites, masked_sites, total_sites = 0, 0, 0
		reMove = False
		num = 0
	
		for record in infile:
			if num == 0:
				for j in range(len(record.seq)):
					check.append(1)
					total_sites = len(check)
			for j in range(len(record.seq)):
				if record.seq[j] not in ["A", "T", "C", "G"]:
					check[j] = 0
				else: 
					''
			unmasked_sites = record.seq.count("A")+record.seq.count("T")+record.seq.count("C")+record.seq.count("G")
			if unmasked_sites < 200:
				reMove = True
			num += 1


		print i
		bases = check.count(1)
		if bases < 150 or reMove:
			print i + " is removed"
		else:
			cmd = "cp "+inDIR+i+" "+outDIR
			os.system(cmd)
		
		handle.close()

	print("--- %s seconds ---" % (time.time() - start_time))	## calculate the running time for the program
		
