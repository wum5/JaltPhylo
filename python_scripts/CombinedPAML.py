import sys
import numpy as np
import pandas as pd
from scipy import stats


def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    

if __name__ == "__main__":
	if len(sys.argv) != 5:
		print "usage: python CombinedPAML.py NSFile PAMLFile FuncFile test(clade/branch)"
		sys.exit()


	outfile = sys.argv[1].split('_')[0]+'_final.txt'
	test = sys.argv[4]

	paml_out = pd.read_table(sys.argv[2], sep='\t')
	paml_out2 = paml_out[paml_out['test.errorstate']=='ok']
	paml_out2 = paml_out2[['gene','alignlength','lrtscore','test.fgdnds2a']] 
	paml_out2.columns = ['gene','size','lrtscore','PAML dN/dS']

	FuncFile = open(sys.argv[3],"r") 
	unames = ['gene', 'function']
	functions = pd.read_table(sys.argv[3], sep='\t', header=None, names=unames)	


	if test == 'branch':
		ns_out = pd.read_table(sys.argv[1], sep='\t') 
		ns_out2 = ns_out[['gene','nonsynonymous_changes','synonymous_changes','ns_ratio']]
		ns_out2.columns = ['gene','NON','SYN','NON/SYN']
		
		data = pd.merge(pd.merge(ns_out2, paml_out2), functions)
		data['pval'] = data['lrtscore'].apply(lambda x : 1 - stats.chi2.cdf(x, 1))
		pList = data['pval'].tolist()
		data['fdr'] = p_adjust_bh(pList)
		data = data.sort_values(by='pval')
		data = data[['gene','size','NON','SYN','NON/SYN','lrtscore','PAML dN/dS','pval','fdr','function']]
		data.to_csv(outfile, sep='\t',index=False)
		

	elif test == 'clade':
		data = pd.merge(paml_out2, functions)
		data['pval'] = data['lrtscore'].apply(lambda x : 1 - stats.chi2.cdf(x, 3))
		pList = data['pval'].tolist()
		data['fdr'] = p_adjust_bh(pList)
		data = data.sort_values(by='PAML dN/dS', ascending=False)
		data = data[['gene','size','lrtscore','PAML dN/dS','pval','fdr','function']]
		data.to_csv(outfile, sep='\t',index=False)
		
		
