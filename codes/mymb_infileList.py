import sys, os


if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: python infileList.py inDIR"
		sys.exit()

	inDIR = sys.argv[1]+"/"
	
	for i in os.listdir(inDIR):
		if i[:4] != "main": continue
		print inDIR+i+"/mymb"
