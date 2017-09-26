#!/usr/bin/env python
"""determine whether the occurrence of the distribution of nonsynonymous
variants is significantly more frequent than the null model by using ms
simulation. Input pattern should be like 00001010101, 1 represents derived 
allele while 0 represents the ancestral alleles"""


import argparse
import numpy as np


def count_nonsyn(mvf_file):
	"""
	count number of nonsynonymous variants from mvf file
	"""
	ns_num = 0
	for line in mvf_file:
		line = line.rstrip()
		if ':' not in line: continue
		codon = line.split(' ')[1]
		if len(codon) > 1:
			ns_num += 1
	return ns_num


def reverse_pattern(pattern):
	rev_pattern = ''
	for x in xrange(len(pattern)):
		if pattern[x] == '0':
			rev_pattern += '1'
		elif pattern[x] == '1':
			rev_pattern += '0'
	return rev_pattern
	

def parse_ms(ms_file, obs_pattern, size):
	"""
	split the simulation dataset into thousands of subset, and 
	calculate the frequency of observed pattern in each subset.
	"""
	pattern = ''
	num, count = 0, 0
	output = []
	
	for line in ms_file:
		if num > int(size):
			output.append(count)
			num, count = 0, 0
		
		line = line.rstrip()
			
		if line == '//' and len(pattern) > 1:
			num += 1
			pattern = pattern[:11]+pattern[12:14]
			if pattern == obs_pattern or reverse_pattern(pattern) == obs_pattern:
				count += 1
			pattern = ''
		elif line == '0' or line == '1':
			pattern += line
			
	return output	


def pval_cal(output, expected_num):
	n = 0
	for x in xrange(len(output)):
		if output[x] >= int(expected_num):
			n += 1 
	pval = float(n+1)/len(output)
	return pval
	

def main():
	parser = argparse.ArgumentParser(prog="ns_count.py", description=__doc__,
	formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--mvf_file', required=True)
	parser.add_argument('-m', '--ms_file', required=True)
	parser.add_argument('-p', '--pattern', required=True)
	parser.add_argument('-n', '--expected_num', required=True)
	args = parser.parse_args()
	
	infile = open(args.mvf_file, "r")
	ns_num = count_nonsyn(infile)	
	print ns_num
	infile.close()
	
	ms_file = open(args.ms_file, "r")
	pattern = args.pattern
	output = parse_ms(ms_file, pattern, ns_num)	
	ms_file.close()

	max_val = max(output)
	min_val = min(output)
	avg_val = np.mean(output)
	
	expected_num = args.expected_num
	pval = pval_cal(output, expected_num)

	print "number of simulations (%d):" % (len(output))
	print "range: %d - %d, mean: %.2f, pval < %f" % (min_val, max_val, avg_val, pval)
	

if __name__ == "__main__":
	main()
	
	
	
