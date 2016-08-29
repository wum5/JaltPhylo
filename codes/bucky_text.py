import os, sys


if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: python addText.py inDir"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	text =  "\nbegin mrbayes;\n\toutgroup Solyc;\n\tset autoclose=yes nowarnings=yes;\n\tlset nst=6 rates=invgamma; [specifies a GTR+I+G]\n\tmcmc ngen=1000000 samplefreq=1000 printfreq=1000 savebrlens=yes filename=mymb;\n\tquit;\nend;"
	
	for i in os.listdir(inDIR):
		with open(inDIR+i, "a") as myfile:
			myfile.write(text)
    		myfile.close()