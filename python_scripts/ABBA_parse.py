import sys, argparse, re
#import matplotlib
#matplotlib.use('Agg')
#from matplotlib import pyplot as plt
#from matplotlib_venn import venn3, venn3_circles
from Bio import SeqIO


def ABBA_parse(file, B2):

	infile = open(file,"r")
	species_num = 0
	species_id = {}
	A1 = 'J.procumbens'
	A2 = 'S.lycopersicum'
	B1 = 'J.calliantha'
	BBAA = []
	ABBA = []
	BABA = []

	for line in infile:
	
		line=line.rstrip()
		if 'J.' in line:
			currID = line.split(' ')[1]
			species_id[currID]=species_num
			species_num+=1
		if A2 in line:
			species_id[A2]=species_num
		if '0:' in line:
			site=line.split(' ')[0]
			nucl=line.split(' ')[1]
			
			if len(nucl)==15:
				x3=nucl[species_id[A1]]
				x4=nucl[species_id[A2]]
				x1=nucl[species_id[B1]]
				x2=nucl[species_id[B2]]
				
				if x1==x2 and x3==x4 and x2!=x3:
					BBAA.append(site)
				elif x1==x3 and x2==x4 and x1!=x2:
					BABA.append(site)
				elif x1==x4 and x2==x3 and x1!=x2:
					ABBA.append(site)
					
			elif '+' in line and nucl[0]!=nucl[1] and nucl[3]==nucl[0]:
				
				if species_id[B2] == 0:
					if re.search('\D+'+str(species_id[A1]),nucl):
						ABBA.append(site)
					elif re.search('\D+'+str(species_id[A2]),nucl):
						BABA.append(site)
					elif re.search('\D+'+str(species_id[B1]),nucl):
						BBAA.append(site)
						
				elif species_id[B1]==0:
					if re.search('\D+'+str(species_id[A1]),nucl):
						BABA.append(site)
					elif re.search('\D+'+str(species_id[A2]),nucl):
						ABBA.append(site)
					elif re.search('\D+'+str(species_id[B2]),nucl):
						BBAA.append(site)

				elif species_id[A1]==0:
					if re.search('\D+'+str(species_id[B1]),nucl):
						BABA.append(site)
					elif re.search('\D+'+str(species_id[B2]),nucl):
						ABBA.append(site)
					elif re.search('\D+'+str(species_id[A2]),nucl):
						BBAA.append(site)

				elif species_id[A2]==0:
					if re.search('\D+'+str(species_id[B2]),nucl):
						BABA.append(site)
					elif re.search('\D+'+str(species_id[B1]),nucl):
						ABBA.append(site)
					elif re.search('\D+'+str(species_id[A1]),nucl):
						BBAA.append(site)


	infile.close()
	return BBAA, BABA, ABBA


def pairwise_dif(abba1, abba2):

	shared=[element for element in abba1 if element in abba2]
	uniq_1=[element for element in abba1 if element not in abba2]
	uniq_2=[element for element in abba2 if element not in abba1]
	return shared, uniq_1, uniq_2


def trio_dif(sp1,sp2,sp3,abba1,abba2,abba3):
	
	YXX = [element for element in abba1 if element not in abba2 and \
	element not in abba3]
	YYX = [element for element in abba1 if element in abba2 and \
	element not in abba3]
	YXY = [element for element in abba1 if element not in abba2 and \
	element in abba3]
	YYY = [element for element in abba1 if element in abba2 and \
	element in abba3]
	XYY = [element for element in abba2 if element in abba3 and \
	element not in abba1]
	XYX = [element for element in abba2 if element not in abba1 and \
	element not in abba3]
	XXY = [element for element in abba3 if element not in abba1 and \
	element not in abba2]
	
	return len(YXX), len(YYX), len(YXY), len(YYY), len(XYY), len(XYX), len(XXY)


def venn_diag(sp1,sp2,sp3,Abc,aBc,ABc,abC,AbC,aBC,ABC):
	
	s = (
    	Abc,   
    	aBc,   
    	ABc,   
    	abC,    
    	AbC,    
    	aBC,  
    	ABC    
	)

	v = venn3(subsets=s, set_labels=(sp1, sp2, sp3))

	v.get_label_by_id('100').set_text(Abc)
	v.get_label_by_id('010').set_text(aBc)
	v.get_label_by_id('110').set_text(ABc)
	v.get_label_by_id('001').set_text(abC)
	v.get_label_by_id('101').set_text(AbC)
	v.get_label_by_id('011').set_text(aBC)
	v.get_label_by_id('111').set_text(ABC)

	v.get_patch_by_id('100').set_color('c')
	v.get_patch_by_id('010').set_color('#993333')
	v.get_patch_by_id('110').set_color('blue')
	
	v.get_patch_by_id('101').set_alpha(0.4)
	v.get_patch_by_id('011').set_alpha(1.0)
	v.get_patch_by_id('111').set_alpha(0.7)
	
	c = venn3_circles(subsets=s, linestyle='solid')
	c[0].set_ls('dotted')  
	c[1].set_ls('dashed')
	c[2].set_lw(1.0)       
	

def trio_plot(infile,sp1_A1,sp2_A1,sp3_A1):

	bbaa1, baba1, abba1 = ABBA_parse(infile, sp1_A1)               
	bbaa2, baba2, abba2 = ABBA_parse(infile, sp2_A1)
	bbaa3, baba3, abba3 = ABBA_parse(infile, sp3_A1)
		
	fig = plt.figure(figsize=(12, 6))
	fig.subplots_adjust(hspace=0.05, wspace=0.05)
	
	ax1 = fig.add_subplot(1, 2, 1)
	YXX,YYX,YXY,YYY,XYY,XYX,XXY = trio_dif(sp1_A1,sp2_A1,sp3_A1,baba1,baba2,baba3)
	venn_diag(sp1_A1,sp2_A1,sp3_A1,YXX,XYX,YYX,XXY,YXY,XYY,YYY)
	plt.title("BABA Pattern")
	ax1.axis('off')
		
	ax2 = fig.add_subplot(1, 2, 2)
	YXX,YYX,YXY,YYY,XYY,XYX,XXY = trio_dif(sp1_A1,sp2_A1,sp3_A1,abba1,abba2,abba3) 
	venn_diag(sp1_A1,sp2_A1,sp3_A1,YXX,XYX,YYX,XXY,YXY,XYY,YYY)
	plt.title("ABBA Pattern")
	ax2.axis('off')

	fig.savefig('%s-%s-%s.pdf' % (sp1_A1,sp2_A1,sp3_A1), bbox_inches='tight')	


def remove_abba_sites(position_list):

	handle = open("/N/dc2/projects/jaltomt/introgression/transcriptome.fa", "r")
	infile = SeqIO.parse(handle, "fasta")
	outfile = open("/N/dc2/projects/jaltomt/introgression/cut_transcriptome.fa","w")
	
	for record in infile:
		myseq = ''
		prev_pos = 0
		for x in xrange(len(position_list)):
			curr_pos = int(position_list[x].split(':')[1])
			myseq += str(record.seq)[prev_pos:curr_pos-1]
			prev_pos = curr_pos
		myseq += str(record.seq)[prev_pos:]
		outfile.write('>'+str(record.id)+'\n')
		outfile.write(myseq+'\n')
			
	infile.close()
	outfile.close()



def main(arguments=sys.argv[1:]):
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-mvf', '--mvf_file', required=True)
	parser.add_argument('-test', '--test_type', choices=["pairwise", "trios", "red-specific"], default="trios",)
	args = parser.parse_args()
	
	infile = open(args.mvf_file	,"r")
	species = ['J.auriculata','J.yungayensis','J.biflora','J.sinuosa','J.aijana',\
	'J.umbellata','J.grandibaccata','J.incahausina','J.dendroidea']
	infile.close()
	
	infile = args.mvf_file
	
	if args.test_type == "trios":
		
		sp1_B2,sp2_B2,sp3_B2='J.umbellata','J.aijana','J.sinuosa'
		trio_plot(infile,sp1_B2,sp2_B2,sp3_B2)
		
		'''		
		sp1_A1,sp2_A1,sp3_A1='J.dendroidea','J.incahausina','J.grandibaccata'
		trio_plot(infile,sp1_A1,sp2_A1,sp3_A1)

		sp1_A1,sp2_A1,sp3_A1='J.yungayensis','J.biflora','J.auriculata'
		trio_plot(infile,sp1_A1,sp2_A1,sp3_A1)
	
		sp1_A1,sp2_A1,sp3_A1='J.biflora','J.sinuosa','J.auriculata'
		trio_plot(infile,sp1_A1,sp2_A1,sp3_A1)
		
		sp1_A1,sp2_A1,sp3_A1='J.grandibaccata','J.dendroidea','J.auriculata'
		trio_plot(infile,sp1_A1,sp2_A1,sp3_A1)
	
		sp1_A1,sp2_A1,sp3_A1='J.grandibaccata','J.incahausina','J.auriculata'
		trio_plot(infile,sp1_A1,sp2_A1,sp3_A1)
	
		sp1_A1,sp2_A1,sp3_A1='J.dendroidea','J.incahausina','J.auriculata'
		trio_plot(infile,sp1_A1,sp2_A1,sp3_A1)
	
		sp1_A1,sp2_A1,sp3_A1='J.aijana','J.umbellata','J.auriculata'
		trio_plot(infile,sp1_A1,sp2_A1,sp3_A1)
		'''
			
	elif args.test_type == 'pairwise':
	
		outfile = open("abba_info.txt", "w")
		outfile.write("species_contrast\tbbaa_common\tbbaa_uniq1\tbbaa_uniq2"
		"\tabba_common\tabba_uniq1\tabba_uniq2\tbaba_common\tbaba_uniq1\t"
		"baba_uniq2\n")

		for x in xrange(len(species)-1):
			sp1_B2=species[x]
			for y in xrange(x+1,len(species)):
				sp2_B2=species[y]
		
				bbaa1, baba1, abba1 = ABBA_parse(infile, sp1_B2)
				bbaa2, baba2, abba2 = ABBA_parse(infile, sp2_B2)

				bbaa_comn, bbaa_uniq1, bbaa_uniq2 = pairwise_dif(bbaa1,bbaa2)
				abba_comn, abba_uniq1, abba_uniq2 = pairwise_dif(abba1,abba2)
				baba_comn, baba_uniq1, baba_uniq2 = pairwise_dif(baba1,baba2)
	
				outfile.write("%s-%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % \
				(sp1_B2,sp2_B2,len(bbaa_comn),len(bbaa_uniq1),len(bbaa_uniq2),len(abba_comn),\
				len(abba_uniq1),len(abba_uniq2),len(baba_comn),len(baba_uniq1),len(baba_uniq2)))
	
		outfile.close()
		
		
	elif args.test_type == "red-specific":
		
		sp1_B2 = 'J.auriculata'
		bbaa1, baba1, abba1 = ABBA_parse(infile, sp1_B2)
		
		for x in xrange(len(species)):
			if species[x] != 'J.auriculata':
				sp2_B2 = species[x]
				bbaa2, baba2, abba2 = ABBA_parse(infile, sp2_B2)

				bbaa_comn, bbaa_uniq1, bbaa_uniq2 = pairwise_dif(bbaa1,bbaa2)
				abba_comn, abba_uniq1, abba_uniq2 = pairwise_dif(abba1,abba2)
				baba_comn, baba_uniq1, baba_uniq2 = pairwise_dif(baba1,baba2)
				
				bbaa1, abba1, baba1 = bbaa_uniq1, abba_uniq1, baba_uniq1
		
		print len(bbaa1), len(abba1), len(baba1)		
		remove_abba_sites(abba1)
	

if __name__=="__main__":
	main()

